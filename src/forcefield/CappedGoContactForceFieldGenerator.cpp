#include "CappedGoContactForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> CappedGoContactForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "cap_val * cap_formula + (1 - cap_val) * lj_formula;"
        "cap_formula = {id}_k * {id}_alpha / roverdist + {id}_k * {id}_beta;"
        "lj_formula = {id}_k * (5 * roverdist^12 - 6 * roverdist^10);"
        "cap_val = step({id}_cap - 1 / roverdist);"
        "roverdist = {id}_r0 / r;",
        fmt::arg("id", ffgen_id_));
    auto contact_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
    contact_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    contact_ff->addPerBondParameter(fmt::format("{}_k", ffgen_id_));
    contact_ff->addPerBondParameter(fmt::format("{}_r0", ffgen_id_));
    contact_ff->addPerBondParameter(fmt::format("{}_alpha", ffgen_id_));
    contact_ff->addPerBondParameter(fmt::format("{}_beta", ffgen_id_));
    contact_ff->addPerBondParameter(fmt::format("{}_cap", ffgen_id_));

    const double cap_slope     = -60 * (std::pow(1 / capping_ratio_, 13.0) - std::pow(1 / capping_ratio_, 11.0));
    const double cap_intercept =  65 * std::pow(1 / capping_ratio_, 12.0) - 66 * std::pow(1 / capping_ratio_, 10.0);

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& indices = indices_vec_[idx];
        contact_ff->addBond(indices.first, indices.second, {ks_[idx], r0s_[idx], cap_slope, cap_intercept, capping_ratio_});
    }

    return contact_ff;
}
