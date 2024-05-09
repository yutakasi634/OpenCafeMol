#include "HarmonicCoMPullingForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> HarmonicCoMPullingForceFieldGenerator::generate() const noexcept
{
    const std::string potential_formula = fmt::format(
        "0.5 * {id}_k * (distance(g1, g2) - {id}_v0)^2",
        fmt::arg("id", ffgen_id_));

    auto com_ff = std::make_unique<OpenMM::CustomCentroidBondForce>(2, potential_formula);

    com_ff->addGlobalParameter(fmt::format("{}_k", ffgen_id_), k_);
    com_ff->addGlobalParameter(fmt::format("{}_v0", ffgen_id_), v0_);
    com_ff->addGroup(first_group_);
    com_ff->addGroup(second_group_);
    com_ff->addBond({0, 1});

    com_ff->setUsesPeriodicBoundaryConditions(use_periodic_);

    return com_ff;
}

