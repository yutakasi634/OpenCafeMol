#ifndef OPEN_AICG2_PLUS_PULLING_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_PULLING_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class PullingForceFieldGenerator : public ForceFieldGeneratorBase
{
  public:
    using parameter_type = std::pair<std::size_t, std::array<double, 3>>;

  public:
    PullingForceFieldGenerator(
        std::vector<parameter_type> params, const bool use_periodic)
        : params_(params), use_periodic_(use_periodic)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
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

    std::string name() const noexcept { return "Pulling"; }

  private:
    std::vector<parameter_type> params_;
    bool                        use_periodic_;
    std::string                 ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_PULLING_FORCE_FIELD_GENERATOR_HPP
