#include "GoContactForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> GoContactForceFieldGenerator::generate() const noexcept
{
    const std::string potential_formula = fmt::format(
        "{id}_k * (5 * ({id}_r0 / r)^12 - 6 * ({id}_r0 / r)^10)",
        fmt::arg("id", ffgen_id_));
    auto contact_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
    contact_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    contact_ff->addPerBondParameter(fmt::format("{}_k", ffgen_id_));
    contact_ff->addPerBondParameter(fmt::format("{}_r0", ffgen_id_));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& indices = indices_vec_[idx];
        contact_ff->addBond(indices.first, indices.second, {ks_[idx], r0s_[idx]});
    }

    return contact_ff;
}
