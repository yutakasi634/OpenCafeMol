#include "GaussianBondForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> GaussianBondForceFieldGenerator::generate() const noexcept
{
    const std::string potential_formula = fmt::format(
        "{id}_k * exp(-(r - {id}_v0)^2 / (2 * {id}_sigma^2))",
        fmt::arg("id", ffgen_id_));

    auto bond_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
    bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    bond_ff->addPerBondParameter(fmt::format("{}_k", ffgen_id_));
    bond_ff->addPerBondParameter(fmt::format("{}_v0", ffgen_id_));
    bond_ff->addPerBondParameter(fmt::format("{}_sigma", ffgen_id_));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& idx_pair = indices_vec_[idx];
        bond_ff->addBond(idx_pair.first, idx_pair.second,
                         {ks_[idx], v0s_[idx], sigmas_[idx]});
    }

    return bond_ff;
}


