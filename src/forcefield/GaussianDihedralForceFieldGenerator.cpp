#include "GaussianDihedralForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> GaussianDihedralForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "{id}_k * exp(-(dt_periodic)^2 /(2*{id}_sigma^2));"
        "dt_periodic = dt - floor((dt + pi)/(2*pi))*(2*pi);"
        "dt = theta-{id}_theta0;"
        "pi = 3.1415926535897932385",
        fmt::arg("id", ffgen_id_));

    auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(potential_formula);
    torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    torsion_ff->addPerTorsionParameter(fmt::format("{}_k", ffgen_id_));
    torsion_ff->addPerTorsionParameter(fmt::format("{}_theta0", ffgen_id_));
    torsion_ff->addPerTorsionParameter(fmt::format("{}_sigma", ffgen_id_));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& indices = indices_vec_[idx];
        torsion_ff->addTorsion(
                indices[0], indices[1], indices[2], indices[3],
                {ks_[idx], theta0s_[idx], sigmas_[idx]});
    }

    return torsion_ff;
}

