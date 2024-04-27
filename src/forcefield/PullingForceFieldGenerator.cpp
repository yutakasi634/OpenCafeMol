#include "PullingForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> PullingForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "-{id}_kx * x - {id}_ky * y - {id}_kz * z",
        fmt::arg("id", ffgen_id_));

    auto pull_ff = std::make_unique<OpenMM::CustomExternalForce>(potential_formula);

    pull_ff->addPerParticleParameter(fmt::format("{}_kx", ffgen_id_));
    pull_ff->addPerParticleParameter(fmt::format("{}_ky", ffgen_id_));
    pull_ff->addPerParticleParameter(fmt::format("{}_kz", ffgen_id_));

    for(const auto& param : params_)
    {
        const std::array<double, 3>& force = param.second;
        pull_ff->addParticle(param.first, {force[0], force[1], force[2]});
    }

    return pull_ff;
}

