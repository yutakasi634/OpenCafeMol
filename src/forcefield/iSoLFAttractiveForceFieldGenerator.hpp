#ifndef OPEN_AICG2_PLUS_ISOLF_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_ISOLF_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP

class iSoLFAttractiveForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    iSoLFAttractiveForceFieldGenerator(
        const std::vector<std::optional<double>> sigmas,
        const std::vector<std::optional<double>> epsilons,
        const std::vector<std::optional<double>> omegas,
        const index_pairs_type& ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : sigmas_(sigmas), epsilons_(epsilons), omegas_(omegas), ignore_list_(ignore_list)
    {
        assert(this->sigmas_.size() == this->epsilons_.size());
        assert(this->sigmas_.size() == this->omegas_.size());

        // make interaction group
        if(ignore_group_pairs.size() == 0)
        {
            std::set<int> participants;
            for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
            {
                if(sigmas_[idx])
                {
                    assert(epsilons_[idx]);
                    assert(omegas_  [idx]);
                    participants.insert(idx);
                }
            }
            interaction_groups_.push_back({ participants, participants });
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
                    assert(epsilons_[idx]);
                    assert(omegas_  [idx]);

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

            interaction_groups_.push_back({ others, others });
            for(const auto& name_group_pair : related_group_map)
            {
                const std::string&   first_name  = name_group_pair.first;
                const std::set<int>& first_group = name_group_pair.second;
                interaction_groups_.push_back({ first_group, others });
                for(const auto& name_group_pair : related_group_map)
                {
                    const std::string&   second_name  = name_group_pair.first;
                    if(!Utility::contains(ignore_group_pairs, { first_name, second_name }))
                    {
                        const std::set<int>& second_group = name_group_pair.second;
                        interaction_groups_.push_back({ first_group, second_group });
                    }
                }
            }
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override
    {
        const std::string potential_formula =
            "-eps *"
            "(step(base-r) +"
            " step(r-base)*step(base+omega-r) * cos(0.5*pi*(r - base)/omega)^2);"
            "base  = sigma*2^(1/6);"
            "sigma = (sigma1+sigma2)*0.5;"
            "eps   = sqrt(eps1*eps2);"
            "omega = (omega1+omega2)*0.5;"
            "pi    = 3.1415926535897932385";
        auto isa_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        isa_ff->addPerParticleParameter("sigma");
        isa_ff->addPerParticleParameter("eps");
        isa_ff->addPerParticleParameter("omega");

        double max_sigma        = std::numeric_limits<double>::min();
        double second_max_sigma = std::numeric_limits<double>::min();
        double max_omega        = std::numeric_limits<double>::min();
        double second_max_omega = std::numeric_limits<double>::min();
        std::set<int> participants;
        for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
        {
            const std::optional<double>& sigma = sigmas_  [idx];
            const std::optional<double>& eps   = epsilons_[idx];
            const std::optional<double>& omega = omegas_  [idx];
            if(sigma && eps && omega)
            {
                const double sigma_val = sigma.value();
                const double omega_val = omega.value();

                isa_ff->addParticle({sigma_val, eps.value(), omega_val});
                participants.insert(idx);

                if(max_sigma <= sigma_val)
                {
                    second_max_sigma = max_sigma;
                    max_sigma        = sigma_val;
                }

                if(max_omega <= omega_val)
                {
                    second_max_omega = max_omega;
                    max_omega        = omega_val;
                }
            }
            else if(!sigma && !eps && !omega)
            {
                isa_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::quiet_NaN()});
            }
            else
            {
                throw std::runtime_error(
                    "[error] iSoLFAttractiveForceFieldGenerator : "
                    "The parameter set is sigma, epsilon and omega. incomplete parameter set "
                    "was given for particle index " + std::to_string(idx) + ".");
            }
        }

        for(const auto& group_pair : interaction_groups_)
        {
            isa_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set cutoff
        isa_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        const double cutoff_distance =
            (max_sigma + second_max_sigma) * 0.5 * std::pow(2.0, 1.0/6.0) +
            (max_omega + second_max_omega) * 0.5;
        isa_ff->setCutoffDistance(cutoff_distance);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            isa_ff->addExclusion(pair.first, pair.second);
        }

        return isa_ff;
    }

  private:
    const std::vector<std::optional<double>> sigmas_;
    const std::vector<std::optional<double>> epsilons_;
    const std::vector<std::optional<double>> omegas_;
    index_pairs_type                         ignore_list_;
    std::vector<interaction_group_type>      interaction_groups_;
};

#endif // OPEN_AICG2_PLUS_ISOLF_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP
