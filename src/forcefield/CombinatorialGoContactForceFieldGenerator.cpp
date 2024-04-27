#include "CombinatorialGoContactForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> CombinatorialGoContactForceFieldGenerator::generate() const
{
    std::string potential_formula = fmt::format(
        "{id}_k *"
            "(5 * ({id}_r0 / distance(a1, d1))^12 -"
            "6 * ({id}_r0 / distance(a1, d1))^10)",
        fmt::arg("id", ffgen_id_));

    auto contact_ff = std::make_unique<OpenMM::CustomHbondForce>(potential_formula);

    contact_ff->addGlobalParameter(fmt::format("{}_k",  ffgen_id_), k_);
    contact_ff->addGlobalParameter(fmt::format("{}_r0", ffgen_id_), r0_);

    // set cutoff
    if(use_periodic_)
    {
        contact_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
    }
    else
    {
        contact_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffNonPeriodic);
    }
    contact_ff->setCutoffDistance(r0_ * cutoff_ratio_);

    for(const auto& donor_particles : donor_indices_vec_)
    {
        contact_ff->addDonor(donor_particles, -1, -1);
    }

    for(const auto& acceptor_particles : acceptor_indices_vec_)
    {
        contact_ff->addAcceptor(acceptor_particles, -1, -1);
    }

    // set exclusion list
    for(const auto& pair : ignore_list_)
    {
        contact_ff->addExclusion(pair.first, pair.second);
    }

    return contact_ff;
}
