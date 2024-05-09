#include "ThreeSPN2ExcludedVolumeForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"

#include <iostream>

ThreeSPN2ExcludedVolumeForceFieldGenerator::ThreeSPN2ExcludedVolumeForceFieldGenerator(
    const double eps, const double cutoff,
    const std::vector<std::optional<double>>& radiuses,
    const index_pairs_type& ignore_list, const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : eps_(eps), cutoff_(cutoff), radiuses_(radiuses), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("TSPN2_EXV{}", ffid.gen()))
{
    interaction_groups_ = extract_interaction_group(radiuses_, ignore_group_pairs, group_vec);
}

std::unique_ptr<OpenMM::Force> ThreeSPN2ExcludedVolumeForceFieldGenerator::generate() const
{
    std::string potential_formula = fmt::format(
        "step(sigma - r)*"
        "(epsilon*((sigma/r)^12 - 2*(sigma/r)^6) + epsilon);"
        "epsilon = sqrt({id}_epsilon1 * {id}_epsilon2);"
        "sigma   = 0.5*({id}_sigma1 + {id}_sigma2);",
        fmt::arg("id", ffgen_id_));

    auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    exv_ff->addPerParticleParameter(fmt::format("{}_epsilon", ffgen_id_));
    exv_ff->addPerParticleParameter(fmt::format("{}_sigma",   ffgen_id_));

    double max_radius        = std::numeric_limits<double>::min();
    double second_max_radius = std::numeric_limits<double>::min();
    for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
    {
        const std::optional<double>& radius = radiuses_[idx];
        if(radius)
        {
            double radius_val = radius.value();
            exv_ff->addParticle({eps_, radius_val});

            if(max_radius < radius_val)
            {
                second_max_radius = max_radius;
                max_radius        = radius_val;
            }
            else if(second_max_radius <= radius_val)
            {
                second_max_radius = radius_val;
            }
        }
        else
        {
            exv_ff->addParticle({
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),});
        }
    }

    // if interaction_groups size is 0, no interaction group will be added,
    // so all the particle in the system will be considerd as participant
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
    std::cerr << "   ThreeSPN2ExcludedVolume : cutoff disntace is "
              << cutoff_distance << " nm" << std::endl;
    exv_ff->setCutoffDistance(cutoff_distance);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        exv_ff->addExclusion(pair.first, pair.second);
    }

    return exv_ff;
}
