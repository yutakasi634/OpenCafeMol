#ifndef OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <limits>
#include <set>
#include <optional>
#include <OpenMM.h>

class ExcludedVolumeForceFieldGenerator final: public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::size_t ffgen_count = 0,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : eps_(eps), cutoff_(cutoff), radiuses_(radiuses), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_count))
    {
        if(ignore_group_pairs.size() == 0)
        {
            // without interaction group force calculation is faster than
            // with that. So if all particle have the parameter for this
            // calculation, force don't use interaction group
            if(std::any_of(radiuses_.begin(), radiuses_.end(),
                          [](const std::optional<double>& val) { return !val; }))
            {
                std::set<int> participants;
                for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
                {
                    if(radiuses_[idx])
                    {
                        participants.insert(idx);
                    }
                }
                interaction_groups_.push_back({ participants, participants });
            }
            else
            {
                std::cerr << "        all particles are participants in this interaction" << std::endl;
            }
        }
        else // make interaction group when group based ignoration specified
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
            for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
            {
                if(radiuses_[idx])
                {
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

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula =
            "EXV"+ffgen_id_str_+"_epsilon*"
            "(sigma_sum)^12*((1/r)^12-EXV"+ffgen_id_str_+"_cutoff_correction);"
            "sigma_sum = EXV"+ffgen_id_str_+"_sigma1 + EXV"+ffgen_id_str_+"_sigma2";
        auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        exv_ff->addPerParticleParameter("EXV"+ffgen_id_str_+"_sigma");
        exv_ff->addGlobalParameter("EXV"+ffgen_id_str_+"_epsilon", eps_);

        double max_radius        = std::numeric_limits<double>::min();
        double second_max_radius = std::numeric_limits<double>::min();
        for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
        {
            const std::optional<double>& radius = radiuses_[idx];
            if(radius)
            {
                exv_ff->addParticle({radius.value()});

                if(max_radius < radius.value())
                {
                    second_max_radius = max_radius;
                    max_radius        = radius.value();
                }
            }
            else
            {
                exv_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
            }
        }

        // if interaction_groups size is 0, no interaction group will be added,
        // so all the particle inthe system will be considerd as participant
        for(const auto& group_pair : interaction_groups_)
        {
            exv_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance   = (max_radius + second_max_radius)*cutoff_;
        const double cutoff_correction = std::pow(1.0 / cutoff_distance, 12);
        std::cerr << "        cutoff disntace is " << cutoff_distance << " nm" << std::endl;
        exv_ff->setCutoffDistance(cutoff_distance);
        exv_ff->addGlobalParameter("EXV"+ffgen_id_str_+"_cutoff_correction", cutoff_correction);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            exv_ff->addExclusion(pair.first, pair.second);
        }

        return exv_ff;
    }

    const std::string name() const noexcept { return "ExcludedVolume"; }

  private:
    double                              eps_;
    double                              cutoff_;
    std::vector<std::optional<double>>  radiuses_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
    const bool                          use_periodic_;
    const std::string                   ffgen_id_str_;
};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
