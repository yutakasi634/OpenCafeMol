#include "iSoLFAttractiveForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"

#include "src/util/Logger.hpp"

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
    log_assert(this->sigmas_.size() == this->epsilons_.size(),
        "iSoLFAttractive: sigma.size(={}) must be the same as epsilon.size(={})",
        sigmas_.size(), epsilons_.size());
    log_assert(this->sigmas_.size() == this->omegas_.size(),
        "iSoLFAttractive: sigma.size(={}) must be the same as omega.size(={})",
        sigmas_.size(), omegas_.size());

    for(std::size_t i=0; i<sigmas_.size(); ++i)
    {
        log_assert(sigmas_.at(i).has_value() == epsilons_.at(i).has_value(),
            "WCA: if particle i(={}) has sigma, it should have epsilon too, and vice versa.", i);
        log_assert(sigmas_.at(i).has_value() == omegas_.at(i).has_value(),
            "WCA: if particle i(={}) has sigma, it should have omega too, and vice versa.", i);
    }

    interaction_groups_ = extract_interaction_group(
            sigmas_, ignore_group_pairs, group_vec);
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

