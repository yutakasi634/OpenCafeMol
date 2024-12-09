#include "DebyeHuckelForceFieldGenerator.hpp"
#include "InteractionGroup.hpp"
#include "ForceFieldIDGenerator.hpp"

#include "src/util/Constants.hpp"

#include <iostream>
#include <iomanip>

DebyeHuckelForceFieldGenerator::DebyeHuckelForceFieldGenerator(
    const double ionic_strength,
    const double temperature, const double cutoff_ratio,
    const std::vector<std::optional<double>>& charges,
    const index_pairs_type& ignore_list,
    const bool use_periodic,
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
    const std::vector<std::optional<std::string>> group_vec)
    : ionic_strength_(ionic_strength), temperature_(temperature),
      cutoff_ratio_(cutoff_ratio), charges_(charges), ignore_list_(ignore_list),
      use_periodic_(use_periodic), ffgen_id_(fmt::format("DH{}", ffid.gen()))
{
    const double epsk = calc_dielectric_water(temperature_, ionic_strength_); // dimensionless
    const double eps0 = Constant::eps0 / Constant::elementary_charge
                                       / Constant::elementary_charge; // [mol/KJ/nm]

    std::cerr << std::scientific << std::setprecision(2)
              << "        eps0 is " << eps0 << " mol/KJ/nm" << std::endl;
    std::cerr << std::fixed
              << "        epsk is " << epsk << " (dimensionless)" << std::endl;

    inv_4_pi_eps0_epsk_ = 1.0 / (4.0 * Constant::pi * eps0 * epsk); // [KJ nm/mol]

    // convert [M] (mol/L) to [mol/nm^3]
    if (!(ionic_strength_ > 0))
    {
        throw std::runtime_error(
                 "[error] ionic strength must be larger than zero");
    }
    
    const double I = ionic_strength_ * 1.0e-24/*[/L]->[/nm^3]*/; // [mol/nm^3]

    debye_length_ =
        std::sqrt((eps0 * epsk * Constant::kB * 1.0e-3 /*[J]->[KJ]*/ * temperature_) / (2. * I)); // [nm]
    std::cerr << "        debye length is " << debye_length_ << " nm" << std::endl;
    abs_cutoff_   = debye_length_ * cutoff_ratio_;
    cutoff_correction_ = std::exp(-cutoff_ratio_) / abs_cutoff_;

    interaction_groups_ = extract_interaction_group(charges_, ignore_group_pairs, group_vec);
}


std::unique_ptr<OpenMM::Force> DebyeHuckelForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "{id}_q1 * {id}_q2 * {id}_inv_4_pi_eps0_epsk * ("
            "exp(-r / {id}_debye_length) / r - {id}_cutoff_correction"
        ")", fmt::arg("id", this->ffgen_id_));

    auto dh_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

    dh_ff->addPerParticleParameter(fmt::format("{}_q", ffgen_id_));
    dh_ff->addGlobalParameter(fmt::format("{}_inv_4_pi_eps0_epsk", ffgen_id_), inv_4_pi_eps0_epsk_); // [KJ nm /mol]
    dh_ff->addGlobalParameter(fmt::format("{}_debye_length", ffgen_id_),       debye_length_); // [nm]
    dh_ff->addGlobalParameter(fmt::format("{}_cutoff_correction", ffgen_id_),  cutoff_correction_);

    for(std::size_t idx=0; idx<charges_.size(); ++idx)
    {
        const std::optional<double>& charge = charges_[idx];
        if(charge)
        {
            dh_ff->addParticle({charge.value()});
        }
        else
        {
            dh_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
        }
    }

    // if interaction_groups size is 0, no interaction group will be added,
    // so all the particle in the system will be considerd as participant
    for(const auto& group_pair : interaction_groups_)
    {
        dh_ff->addInteractionGroup(group_pair.first, group_pair.second);
    }

    // set pbc condition
    if(use_periodic_)
    {
        dh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
    }
    else
    {
        dh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
    }

    // set cutoff
    std::cerr << "    DebyeHuckel                   : cutoff distance is "
              << abs_cutoff_ << " nm" << std::endl;
    dh_ff->setCutoffDistance(abs_cutoff_);

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        dh_ff->addExclusion(pair.first, pair.second);
    }

    return dh_ff;
}

