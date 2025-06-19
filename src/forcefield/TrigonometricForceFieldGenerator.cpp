#include "TrigonometricForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"

#include "src/util/Logger.hpp"

#include <iostream>

TrigonometricForceFieldGenerator::TrigonometricForceFieldGenerator(
    const std::vector<std::optional<double>>& epsilons,
    const std::vector<std::optional<double>>& sigmas,
    const std::vector<std::optional<double>>& omegas,
    const index_pairs_type& ignore_list, const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : epsilons_(epsilons), sigmas_(sigmas), omegas_(omegas),
      ignore_list_(ignore_list), use_periodic_(use_periodic),
      ffgen_id_(fmt::format(fmt::format("TRI{}", ffid.gen())))
{
    log_assert(this->epsilons_.size() == this->sigmas_.size(),
            "Trigonometric: sigmas_.size(={}) must be the same as"
            " epsilons_.size(={})", sigmas_.size(), epsilons_.size());
    log_assert(this->epsilons_.size() == this->omegas_.size(),
            "Trigonometric: omegas_.size(={}) must be the same as"
            " epsilons_.size(={})", omegas_.size(), epsilons_.size());

    for(std::size_t i=0; i<sigmas_.size(); ++i)
    {
        log_assert(epsilons_.at(i).has_value() == sigmas_.at(i).has_value(),
                "Trigonometric: if particle i(={}) has epsilon, "
                "it should have sigma too, and vice versa.", i);
        log_assert(epsilons_.at(i).has_value() == omegas_.at(i).has_value(),
                "Trigonometric: if particle i(={}) has epsilon, "
                "it should have omega too, and vice versa.", i);
    }

    interaction_groups_ = extract_interaction_group(
            this->epsilons_, ignore_group_pairs, group_vec);
}

std::unique_ptr<OpenMM::Force> TrigonometricForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "-epsilon *"
        "    (step(threthold+omega-r) +"
        "      step(r-threthold)*step(threthold+omega-r) *"
        "          (2*r_th_omeg^3 - 3*r_th_omeg^2));"
        "r_th_omeg = (r-threthold)/omega;"
        "threthold = sigma*2^(1/6);"
        "epsilon   = ({id}_epsilon1 * {id}_epsilon2)^(1/2);"
        "sigma     = ({id}_sigma1 + {id}_sigma2) * 0.5;"
        "omega     = ({id}_omega1 + {id}_omega2) * 0.5",
        fmt::arg("id", ffgen_id_));
    auto tri_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    tri_ff->addPerParticleParameter(fmt::format("{}_epsilon", ffgen_id_));
    tri_ff->addPerParticleParameter(fmt::format("{}_sigma", ffgen_id_));
    tri_ff->addPerParticleParameter(fmt::format("{}_omega", ffgen_id_));

    double max_radius        = std::numeric_limits<double>::min();
    double second_max_radius = std::numeric_limits<double>::min();
    for(std::size_t idx=0; idx<sigmas_.size(); ++idx)
    {
        const std::optional<double>& epsilon = epsilons_[idx];
        const std::optional<double>& sigma   = sigmas_[idx];
        const std::optional<double>& omega   = omegas_[idx];
        if(epsilon && sigma && omega)
        {
            double sigma_val = sigma.value();
            double omega_val = omega.value();
            const double radius = sigma_val + omega_val;
            tri_ff->addParticle(
                    {epsilon.value(), sigma_val, omega_val});

            if(max_radius <= radius)
            {
                second_max_radius = max_radius;
                max_radius        = radius;
            }
            else if(second_max_radius <= radius)
            {
                second_max_radius = radius;
            }
        }
        else if(!epsilon && !sigma && !omega)
        {
            tri_ff->addParticle({std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN()});
        }
        else
        {
            throw std::runtime_error(
                "[error] TrigonometricForceFieldGenerator : "
                "The parameter set is epsilon, sigma and omega. incomplete parameter set "
                "was given for particle idx " + std::to_string(idx) + ".");
        }
    }

    // if interaction_group size is 0, no interaction group will be added,
    // so all the particle in the system will be considered as participant
    for(const auto& group_pair : interaction_groups_)
    {
        tri_ff->addInteractionGroup(group_pair.first, group_pair.second);
    }


    // set pbc condition
    if(use_periodic_)
    {
        tri_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
    }
    else
    {
        tri_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    // set cutoff
    const double cutoff_distance =
        (max_radius + second_max_radius) * 0.5;
    std::cerr << "    Trigonomeric                  : cutoff distance is "
              << cutoff_distance << " nm" << std::endl;
    tri_ff->setCutoffDistance(cutoff_distance);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        tri_ff->addExclusion(pair.first, pair.second);
    }

    return tri_ff;
}
