#include "PositionRestraintForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <iostream>

std::unique_ptr<OpenMM::Force> PositionRestraintForceFieldGenerator::generate() const
{
    std::string potential_formula = fmt::format(
            "{id}_k * (((x-{id}_x0)^2 + (y-{id}_y0)^2 + (z-{id}_z0)^2)^(1/2) - {id}_v0)^2",
            fmt::arg("id", ffgen_id_));

    if(use_periodic_)
    {
        std::string potential_formula = fmt::format(
                "{id}_k * (periodicdistance(x, y, z, {id}_x, {id}_y, {id}_z) - {id}_v0)^2",
                fmt::arg("id", ffgen_id_));
    }

    auto pr_ff = std::make_unique<OpenMM::CustomExternalForce>(potential_formula);

    pr_ff->addPerParticleParameter(fmt::format("{}_x0", ffgen_id_));
    pr_ff->addPerParticleParameter(fmt::format("{}_y0", ffgen_id_));
    pr_ff->addPerParticleParameter(fmt::format("{}_z0", ffgen_id_));
    pr_ff->addPerParticleParameter(fmt::format("{}_k",  ffgen_id_));
    pr_ff->addPerParticleParameter(fmt::format("{}_v0", ffgen_id_));

    for(std::size_t param_idx=0; param_idx<indices_.size(); ++param_idx)
    {
        const auto   pos = positions_[param_idx];
        const double k   = ks_[param_idx];
        const double v0  = v0s_[param_idx];
        pr_ff->addParticle(indices_[param_idx], {pos[0], pos[1], pos[2], k, v0});
    }

    return pr_ff;
}

