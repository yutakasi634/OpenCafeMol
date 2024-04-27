#include "iSoLFAttractiveForceFieldGenerator.hpp"

#include "ForceFieldIDGenerator.hpp"
#include "src/util/Utility.hpp"

#include <algorithm>
#include <iostream>

iSoLFAttractiveForceFieldGenerator::iSoLFAttractiveForceFieldGenerator(
    const std::vector<std::optional<double>> sigmas,
    const std::vector<std::optional<double>> epsilons,
    const std::vector<std::optional<double>> omegas,
    const index_pairs_type& ignore_list,
    const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : sigmas_(sigmas), epsilons_(epsilons), omegas_(omegas), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("iSA{}", ffid.gen()))
{
    assert(this->sigmas_.size() == this->epsilons_.size());
    assert(this->sigmas_.size() == this->omegas_.size());

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
                    assert(omegas_  [idx]);
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

std::unique_ptr<OpenMM::Force> iSoLFAttractiveForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "-eps *"
        "(step(base-r) +"
        " step(r-base)*step(base+omega-r) * cos(0.5*pi*(r - base)/omega)^2);"
        "base  = sigma * 2^(1/6);"
        "sigma = ({id}_sigma1 + {id}_sigma2) * 0.5;"
        "eps   = sqrt({id}_eps1 * {id}_eps2);"
        "omega = ({id}_omega1 + {id}_omega2) * 0.5;"
        "pi    = 3.1415926535897932385",
        fmt::arg("id", ffgen_id_));

    auto isa_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    isa_ff->addPerParticleParameter(fmt::format("{}_sigma", ffgen_id_));
    isa_ff->addPerParticleParameter(fmt::format("{}_eps",   ffgen_id_));
    isa_ff->addPerParticleParameter(fmt::format("{}_omega", ffgen_id_));

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
            else if(second_max_sigma <= sigma_val)
            {
                second_max_sigma = sigma_val;
            }

            if(max_omega <= omega_val)
            {
                second_max_omega = max_omega;
                max_omega        = omega_val;
            }
            else if(second_max_omega <= omega_val)
            {
                second_max_omega = omega_val;
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

    // if interaction_groups size is 0, no interaction group will be added,
    // so all the particle in the system will be considerd as participant
    for(const auto& group_pair : interaction_groups_)
    {
        isa_ff->addInteractionGroup(group_pair.first, group_pair.second);
    }

    // set pbc condition
    if(use_periodic_)
    {
        isa_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
    }
    else
    {
        isa_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    // set cutoff
    const double cutoff_distance =
        (max_sigma + second_max_sigma) * 0.5 * std::pow(2.0, 1.0/6.0) +
        (max_omega + second_max_omega) * 0.5;
    std::cerr << "    iSoLFAttractive               : cutoff disntace is "
              << cutoff_distance << " nm" << std::endl;
    isa_ff->setCutoffDistance(cutoff_distance);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        isa_ff->addExclusion(pair.first, pair.second);
    }

    return isa_ff;
}

