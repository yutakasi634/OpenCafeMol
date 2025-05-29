#include "ExcludedVolumeForceFieldGenerator.hpp"
#include "InteractionGroup.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <iostream>
#include <limits>

// The size of the vector representing the per-particle parameters
// (in this case, radii) must match the system size. This must be guaranteed
// inside the read function.
ExcludedVolumeForceFieldGenerator::ExcludedVolumeForceFieldGenerator(
    const double eps, const double cutoff,
    const std::vector<std::optional<double>>& radii,
    const index_pairs_type& ignore_list, const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : eps_(eps), cutoff_(cutoff), radii_(radii), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("EXV{}", ffid.gen()))
{
    interaction_groups_ = extract_interaction_group(radii_, ignore_group_pairs, group_vec);
}

std::unique_ptr<OpenMM::Force> ExcludedVolumeForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "{id}_epsilon * (sigma_sum)^12 * ((1/r)^12 - {id}_cutoff_correction);"
        "sigma_sum = {id}_sigma1 + {id}_sigma2", fmt::arg("id", ffgen_id_));

    auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    exv_ff->addPerParticleParameter(fmt::format("{}_sigma", ffgen_id_));
    exv_ff->addGlobalParameter(fmt::format("{}_epsilon", ffgen_id_), eps_);

    double max_radius        = std::numeric_limits<double>::min();
    double second_max_radius = std::numeric_limits<double>::min();
    for(std::size_t idx=0; idx<radii_.size(); ++idx)
    {
        const std::optional<double>& radius = radii_[idx];
        if(radius)
        {
            const double radius_val = radius.value();
            exv_ff->addParticle({radius_val});

            if(max_radius < radius_val)
            {
                second_max_radius = max_radius;
                max_radius        = radius_val;
            }
            else if(second_max_radius < radius_val)
            {
                second_max_radius = radius_val;
            }
        }
        else
        {
            exv_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
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
    const double cutoff_correction = std::pow(1.0 / cutoff_distance, 12);
    std::cerr << "    ExcludedVolume                : cutoff disntace is "
              << cutoff_distance << " nm" << std::endl;
    exv_ff->setCutoffDistance(cutoff_distance);
    exv_ff->addGlobalParameter(fmt::format("{}_cutoff_correction", ffgen_id_), cutoff_correction);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        exv_ff->addExclusion(pair.first, pair.second);
    }

    return exv_ff;
}


