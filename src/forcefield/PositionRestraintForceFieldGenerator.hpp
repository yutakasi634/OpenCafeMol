#ifndef OPEN_AICG2_PLUS_POSITION_RESTRAINT_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_POSITION_RESTRAINT_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class PositionRestraintForceFieldGenerator : public ForceFieldGeneratorBase
{
  public:
    PositionRestraintForceFieldGenerator(
        std::vector<std::size_t> indices, std::vector<std::array<double, 3>> positions,
        std::vector<double>      ks,      std::vector<double>                v0s)
        : indices_(indices), positions_(positions), ks_(ks), v0s_(v0s),
          ffgen_id_(fmt::format("PR{}", ffid.gen()))
    {
        assert(this->indices_.size() == this->positions_.size());
        assert(this->indices_.size() == this->ks_.size());
        assert(this->indices_.size() == this->v0s_.size());
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula = fmt::format(
                "{id}_k * (periodicdistance(x, y, z, {id}_x0, {id}_y0, {id}_z0) - {id}_v0)^2",
                fmt::arg("id", ffgen_id_));

        auto pr_ff = std::make_unique<OpenMM::CustomExternalForce>(potential_formula);

        pr_ff->addPerParticleParameter(fmt::format("{}_x0", ffgen_id_));
        pr_ff->addPerParticleParameter(fmt::format("{}_y0", ffgen_id_));
        pr_ff->addPerParticleParameter(fmt::format("{}_z0", ffgen_id_));
        pr_ff->addPerParticleParameter(fmt::format("{}_k",  ffgen_id_));
        pr_ff->addPerParticleParameter(fmt::format("{}_v0", ffgen_id_));

        for(std::size_t param_idx=0; param_idx<indices_.size(); param_idx++)
        {
            const auto   pos = positions_[param_idx];
            const double k   = ks_[param_idx];
            const double v0  = v0s_[param_idx];
            pr_ff->addParticle(indices_[param_idx], {pos[0], pos[1], pos[2], k, v0});
        }

        return pr_ff;
    }

    std::string name() const noexcept override { return "PositionRestraint"; }

  private:
    std::vector<std::size_t>           indices_;
    std::vector<std::array<double, 3>> positions_;
    std::vector<double>                ks_;
    std::vector<double>                v0s_;
    std::string                        ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_POSITION_RESTRAINT_FORCE_FIELD_GENERATOR_HPP
