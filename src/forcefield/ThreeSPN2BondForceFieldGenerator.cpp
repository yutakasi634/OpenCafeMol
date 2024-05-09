#include "ThreeSPN2BondForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> ThreeSPN2BondForceFieldGenerator::generate() const noexcept
{
    std::string potential_formula = fmt::format(
        "{id}_k2 * (r - {id}_v0)^2 + 100.0 * {id}_k4 * (r - {id}_v0)^4",
        fmt::arg("id", ffgen_id_));

    auto bond_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);

    bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    bond_ff->addPerBondParameter(fmt::format("{}_k2", ffgen_id_));
    bond_ff->addPerBondParameter(fmt::format("{}_k4", ffgen_id_));
    bond_ff->addPerBondParameter(fmt::format("{}_v0", ffgen_id_));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& idx_pair = indices_vec_[idx];
        bond_ff->addBond(idx_pair.first, idx_pair.second,
                         {k2s_[idx], k4s_[idx], v0s_[idx]});
    }

    return bond_ff;
}

