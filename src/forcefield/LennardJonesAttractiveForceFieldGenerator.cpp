#include "LennardJonesAttractiveForceFieldGenerator.hpp"

#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"
#include "src/util/Logger.hpp"

#include <iostream>

LennardJonesAttractiveForceFieldGenerator::LennardJonesAttractiveForceFieldGenerator(const double cutoff_ratio,
    const std::vector<std::optional<double>> epsilons,
    const std::vector<std::optional<double>> sigmas,
    const index_pairs_type& ignore_list, const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : cutoff_ratio_(cutoff_ratio), epsilons_(epsilons),
      sigmas_(sigmas), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("LJAT{}", ffid.gen()))
{
    assert(this->epsilons_.size() == this->sigmas_.size());

    log_assert(this->sigmas_.size() == this->epsilons_.size(),
        "LennardJonesAttractive: sigma.size(={}) must be the same as epsilon.size(={})",
        sigmas_.size(), epsilons_.size());

    for(std::size_t i=0; i<sigmas_.size(); ++i)
    {
        log_assert(sigmas_.at(i).has_value() == epsilons_.at(i).has_value(),
            "LennardJonesAttractive: if particle i(={}) has sigma, "
            "it should have epsilon too, and vice versa.", i);
    }

    interaction_groups_ = extract_interaction_group(
            this->epsilons_, ignore_group_pairs, group_vec);

    cutoff_correction_ =
        4.0*(std::pow(1.0 / cutoff_ratio, 12) - std::pow(1.0 / cutoff_ratio, 6));
}

std::unique_ptr<OpenMM::Force> LennardJonesAttractiveForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "epsilon *"
        "    (step(r-threthold) * 4 * (sigma_r_12 - sigma_r_6) -"
        "     step(threthold-r));"
        "sigma_r_12 = sigma_r_6^2;"
        "sigma_r_6  = sigma_r^6;"
        "sigma_r    = sigma/r;"
        "threthold  = sigma*2^(1/6);"
        "sigma      = ({id}_sigma1 + {id}_sigma2) * 0.5;"
        "epsilon    = ({id}_epsilon1 + {id}_epsilon2) * 0.5",
        fmt::arg("id", ffgen_id_));
    auto ljattr_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    ljattr_ff->addPerParticleParameter(fmt::format("{}_epsilon", ffgen_id_));
    ljattr_ff->addPerParticleParameter(fmt::format("{}_sigma",   ffgen_id_));

    double max_sigma        = std::numeric_limits<double>::min();
    double second_max_sigma = std::numeric_limits<double>::min();
    for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
    {
        const std::optional<double>& epsilon = epsilons_[idx];
        const std::optional<double>& sigma   = sigmas_[idx];
        if(sigma && epsilon)
        {
            double sigma_val = sigma.value();
            ljattr_ff->addParticle({epsilon.value(), sigma_val});

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
            ljattr_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN()});
        }
        else
        {
            throw std::runtime_error(
                "[error] LennardJonesAttractiveForceFieldGenerator : "
                "The parameter set is epsilon and sigma. incomplete parameter set "
                "was given for particle idx " + std::to_string(idx) + ".");
        }
    }

        // if interaction_group size is 0, no interaction group will be added,
        // so all the particle in the system will be considered as participant
        for(const auto& group_pair : interaction_groups_)
        {
            ljattr_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            ljattr_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            ljattr_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance =
            (max_sigma + second_max_sigma) * 0.5 * cutoff_ratio_;
        std::cerr << "    LennardJonesAttractive        : cutoff distance is "
                  << cutoff_distance << " nm" << std::endl;
        ljattr_ff->setCutoffDistance(cutoff_distance);
        ljattr_ff->setUseSwitchingFunction(true);
        ljattr_ff->setSwitchingDistance(0.9*cutoff_distance);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            ljattr_ff->addExclusion(pair.first, pair.second);
        }

        return ljattr_ff;
    }

