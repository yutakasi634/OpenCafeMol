#include "CosineDihedralForceFieldGenerator.hpp"

#include <regex>

std::unique_ptr<OpenMM::Force> CosineDihedralForceFieldGenerator::generate() const noexcept
{
    std::string potential_formula =  "k_periodic * (1 - cos(n0 * (theta - t0)))";
    // The 3SPN2C DNA model paper (Freeman et al., JCP, 2014) shows k*(1 + cos x) formula.
    // However, we should note that this is a misprint of k*(1 - cos x).
    // Indeed, the other implementations of 3SPN2 such as open3spn2, LAMMPS, and CafeMol
    // use k* (1 - cos x) formula.

    const std::map<std::string, std::string> ff_params =
    {
        {"k_periodic", ffgen_id_ + "_k_periodic"},
        {"t0",         ffgen_id_ + "_t0"},
        {"n0",         ffgen_id_ + "_n0"},
    };

    for(const auto& param : ff_params)
    {
          potential_formula = std::regex_replace(
            potential_formula, std::regex(param.first), param.second);
    }

    auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(potential_formula);

    torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);

    torsion_ff->addPerTorsionParameter(ff_params.at("k_periodic"));
    torsion_ff->addPerTorsionParameter(ff_params.at("t0"));
    torsion_ff->addPerTorsionParameter(ff_params.at("n0"));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& indices = indices_vec_[idx];
        torsion_ff->addTorsion(
                indices[0], indices[1], indices[2], indices[3],
                {ks_[idx], theta0s_[idx], ns_[idx]});
    }

    return torsion_ff;
}
