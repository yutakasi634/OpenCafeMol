#include "WeeksChandlerAndersenForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "InteractionGroup.hpp"

#include "src/util/Logger.hpp"

#include <iostream>

WeeksChandlerAndersenForceFieldGenerator::WeeksChandlerAndersenForceFieldGenerator(
    const std::vector<std::optional<double>> sigmas,
    const std::vector<std::optional<double>> epsilons,
    const index_pairs_type& ignore_list,
    const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : sigmas_(sigmas), epsilons_(epsilons), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("WCA{}", ffid.gen()))
{
    log_assert(this->sigmas_.size() == this->epsilons_.size(),
        "WCA: sigma.size(={}) must be the same as epsilon.size(={})",
        sigmas_.size(), epsilons_.size());

    for(std::size_t i=0; i<sigmas_.size(); ++i)
    {
        log_assert(sigmas_.at(i).has_value() == epsilons_.at(i).has_value(),
            "WCA: if particle i(={}) has sigma, it should have epsilon too, and vice versa.", i);
    }
    interaction_groups_ = extract_interaction_group(
            sigmas_, ignore_group_pairs, group_vec);
}

std::unique_ptr<OpenMM::Force> WeeksChandlerAndersenForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "step(sigma*2^(1/6)-r) * (4*eps * (sigma_r_12 - sigma_r_6) + eps);"
        "sigma_r_12 = sigma_r_6^2;"
        "sigma_r_6 = sigma_r^6;"
        "sigma_r = sigma/r;"
        "eps     = sqrt({id}_eps1 * {id}_eps2);"
        "sigma   = ({id}_sigma1 + {id}_sigma2) * 0.5",
        fmt::arg("id", ffgen_id_));

    auto wca_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    wca_ff->addPerParticleParameter(fmt::format("{}_sigma", ffgen_id_));
    wca_ff->addPerParticleParameter(fmt::format("{}_eps", ffgen_id_));

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
            else if(second_max_sigma <= sigma_val)
            {
                second_max_sigma = sigma_val;
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

    // if interaction_group size is 0, no interaction group will be added,
    // so all the particle in the system will be considerd as participant
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
    std::cerr << "    WeeksChandlerAndersen         : cutoff disntace is "
              << cutoff_distance << " nm" << std::endl;
    wca_ff->setCutoffDistance(cutoff_distance);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        wca_ff->addExclusion(pair.first, pair.second);
    }

    return wca_ff;
}

