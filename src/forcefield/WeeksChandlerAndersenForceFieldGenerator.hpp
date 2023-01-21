#ifndef OPEN_AICG2_PLUS_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>

class WeeksChandlerAndersenForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    WeeksChandlerAndersenForceFieldGenerator(
        const std::vector<std::optional<double>> sigmas,
        const std::vector<std::optional<double>> epsilons,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : sigmas_(sigmas), epsilons_(epsilons), ignore_list_(ignore_list),
          use_periodic_(use_periodic)
    {
        assert(this->sigmas_.size() == this->epsilons_.size());

        // make interaction group
        if(ignore_group_pairs.size() == 0)
        {
            // without interaction group force calculation is faster than
            // with that. So if all particle have the parameter for this
            // calculation, force don't use interaction group
            if(std::any_of(sigmas_.begin(), sigmas_.end(),
                          [](const std::optional<double>& val) { return !val; }))
            {
                std::set<int> participants;
                for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
                {
                    if(sigmas_[idx])
                    {
                        assert(epsilons_[idx]);
                        participants.insert(idx);
                    }
                }
                interaction_groups_.push_back({ participants, participants });
            }
        }
        else // group based ignoration specified case
        {
            std::set<std::string> related_group_names;
            for(const auto& name_group_pair : ignore_group_pairs)
            {
                related_group_names.insert(name_group_pair.first);
                related_group_names.insert(name_group_pair.second);
            }

            std::map<std::string, std::set<int>> related_group_map;
            for(const auto& name : related_group_names)
            {
                related_group_map.insert(std::make_pair(name, std::set<int>()));
            }

            std::set<int> others;
            for(std::size_t idx=0; idx<sigmas.size(); ++idx)
            {
                if(sigmas_[idx])
                {
                    assert(epsilons[idx]);

                    if(group_vec[idx])
                    {
                        const std::string group_name = group_vec[idx].value();
                        if(related_group_names.count(group_name) != 0)
                        {
                            related_group_map.at(group_name).insert(idx);
                        }
                        else
                        {
                            others.insert(idx);
                        }
                    }
                    else
                    {
                        others.insert(idx);
                    }
                }
            }

            std::vector<std::pair<std::string, std::set<int>>> related_group_vec;
            for(const auto& name_group_pair : related_group_map)
            {
                related_group_vec.push_back(name_group_pair);
            }

            interaction_groups_.push_back({ others, others });
            for(std::size_t idx_i=0 ; idx_i<related_group_vec.size(); ++idx_i)
            {
                const auto& name_group_pair_i = related_group_vec[idx_i];
                const std::string&   first_name  = name_group_pair_i.first;
                const std::set<int>& first_group = name_group_pair_i.second;
                interaction_groups_.push_back({ first_group, others });
                for(std::size_t idx_j=idx_i; idx_j<related_group_vec.size(); ++idx_j)
                {
                    const auto& name_group_pair_j = related_group_vec[idx_j];
                    const std::string&   second_name  = name_group_pair_j.first;
                    if(!Utility::contains(ignore_group_pairs, { first_name, second_name }))
                    {
                        const std::set<int>& second_group = name_group_pair_j.second;
                        interaction_groups_.push_back({ first_group, second_group });
                    }
                }
            }
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override
    {
        const std::string potential_formula =
            "step(sigma*2^(1/6)-r) * (4*eps * (sigma_r_12 - sigma_r_6) + eps);"
            "sigma_r_12 = sigma_r_6^2;"
            "sigma_r_6 = sigma_r^6;"
            "sigma_r = sigma/r;"
            "eps     = sqrt(eps1*eps2);"
            "sigma   = (sigma1+sigma2)*0.5";
        auto wca_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        wca_ff->addPerParticleParameter("sigma");
        wca_ff->addPerParticleParameter("eps");

        double max_sigma        = std::numeric_limits<double>::min();
        double second_max_sigma = std::numeric_limits<double>::min();
        for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
        {
            const std::optional<double>& sigma   = sigmas_[idx];
            const std::optional<double>& epsilon = epsilons_[idx];
            if(sigma && epsilon)
            {
                double sigma_val = sigma.value();
                wca_ff->addParticle({sigma_val, epsilon.value()});

                if(max_sigma <= sigma_val)
                {
                    second_max_sigma = max_sigma;
                    max_sigma        = sigma_val;
                }
            }
            else if(!sigma && !epsilon)
            {
                wca_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()});
            }
            else
            {
                throw std::runtime_error(
                    "[error] WeeksChandlerAndersenForceFieldGenerator : "
                    "The parameter set is sigma and epsilon. incomplete parameter set "
                    "was given for particle idx " + std::to_string(idx) + ".");
            }
        }

        // if interaction_groups size is 0, no interaction group will be added,
        // so all the particle inthe system will be considerd as participant
        for(const auto& group_pair : interaction_groups_)
        {
            wca_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            wca_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            wca_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance =
            (max_sigma + second_max_sigma) * 0.5 * std::pow(2.0, 1.0/6.0);
        std::cerr << "        cutoff disntace is " << cutoff_distance << " nm" << std::endl;
        wca_ff->setCutoffDistance(cutoff_distance);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            wca_ff->addExclusion(pair.first, pair.second);
        }

        return wca_ff;
    }

  private:
    const std::vector<std::optional<double>> sigmas_;
    const std::vector<std::optional<double>> epsilons_;
    index_pairs_type                         ignore_list_;
    std::vector<interaction_group_type>      interaction_groups_;
    const bool                               use_periodic_;
};

#endif // OPEN_AICG2_PLUS_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
