#include "HarmonicBondForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> HarmonicBondForceFieldGenerator::generate() const noexcept
{
    auto bond_ff = std::make_unique<OpenMM::HarmonicBondForce>();
    bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& idx_pair = indices_vec_[idx];
        bond_ff->addBond(idx_pair.first, idx_pair.second, v0s_[idx], ks_[idx]);
    }

    return bond_ff;
}
