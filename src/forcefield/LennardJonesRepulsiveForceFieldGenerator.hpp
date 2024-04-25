#ifndef OPEN_AICG2_PLUS_LENNARD_JONES_REPULSIVE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_LENNARD_JONES_REPULSIVE_FORCE_FIELD_GENERATOR_HPP

#include "src/util/Utility.hpp"
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

class LennardJonesRepulsiveForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    LennardJonesRepulsiveForceFieldGenerator(const double cutoff_ratio,
        const std::vector<std::optional<double>> epsilons,
        const std::vector<std::optional<double>> sigmas,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : cutoff_ratio_(cutoff_ratio), epsilons_(epsilons),
          sigmas_(sigmas), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("LJRP{}", ffid.gen()))
    {
        assert(this->epsilons_.size() == this->sigmas_.size());

        // make interaction group
        if(ignore_group_pairs.size() == 0)
        {
            // without interaction group force calculation is faster than
            // with that. So if all particle have the parameter for this
            // calculation, force don't use interaction group
            if(std::any_of(epsilons_.begin(), epsilons_.end(),
                          [](const std::optional<double>& val) { return !val; }))
            {
                std::set<int> participants;
                for(std::size_t idx=0; idx<epsilons_.size(); ++idx)
                {
                    if(epsilons_[idx])
                    {
                        assert(sigmas_[idx]);
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
            for(std::size_t idx=0; idx<epsilons_.size(); ++idx)
            {
                if(epsilons_[idx])
                {
                    assert(sigmas_[idx]);

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
            for(std::size_t idx_i=0; idx_i<related_group_vec.size(); ++idx_i)
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
        const std::string potential_formula = fmt::format(
            "step(threthold-r) * epsilon * (4*(sigma_r_12 - sigma_r_6) + 1);"
            "sigma_r_12 = sigma_r_6^2;"
            "sigma_r_6  = sigma_r^6;"
            "sigma_r    = sigma/r;"
            "threthold  = sigma*2^(1/6);"
            "sigma      = ({id}_sigma1 + {id}_sigma2) * 0.5;"
            "epsilon    = ({id}_epsilon1 + {id}_epsilon2) * 0.5",
            fmt::arg("id", ffgen_id_));
        auto ljrepu_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        ljrepu_ff->addPerParticleParameter(fmt::format("{}_epsilon", ffgen_id_));
        ljrepu_ff->addPerParticleParameter(fmt::format("{}_sigma",   ffgen_id_));

        double max_sigma        = std::numeric_limits<double>::min();
        double second_max_sigma = std::numeric_limits<double>::min();
        for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
        {
            const std::optional<double>& epsilon = epsilons_[idx];
            const std::optional<double>& sigma   = sigmas_[idx];
            if(sigma && epsilon)
            {
                double sigma_val = sigma.value();
                ljrepu_ff->addParticle({epsilon.value(), sigma_val});

                if(max_sigma <= sigma_val)
                {
                    second_max_sigma = max_sigma;
                    max_sigma        = sigma_val;
                }
                else if(second_max_sigma <= sigma_val)
                {
                    second_max_sigma = sigma_val;
                }
            }
            else if(!sigma && !epsilon)
            {
                ljrepu_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN()});
            }
            else
            {
                throw std::runtime_error(
                    "[error] LennardJonesRepulsiveForceFieldGenerator : "
                    "The parameter set is epsilon and sigma. incomplete parameter set "
                    "was given for particle idx " + std::to_string(idx) + ".");
            }
        }

        // if interaction_group size is 0, no interaction group will be added,
        // so all the particle in the system will be considered as participant
        for(const auto& group_pair : interaction_groups_)
        {
            ljrepu_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            ljrepu_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            ljrepu_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance =
            (max_sigma + second_max_sigma) * 0.5 * cutoff_ratio_;
        std::cerr << "    LennardJonesRepulsive         : cutoff distance is "
                  << cutoff_distance << " nm" << std::endl;
        ljrepu_ff->setCutoffDistance(cutoff_distance);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            ljrepu_ff->addExclusion(pair.first, pair.second);
        }

        return ljrepu_ff;
    }

    std::string name() const noexcept { return "LennardJonesRepulsive"; }

  private:
    double                              cutoff_ratio_;
    std::vector<std::optional<double>>  epsilons_;
    std::vector<std::optional<double>>  sigmas_;
    index_pairs_type                    ignore_list_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;
    std::vector<interaction_group_type> interaction_groups_;
};


#endif // OPEN_AICG2_PLUS_LENNARD_JONES_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP
