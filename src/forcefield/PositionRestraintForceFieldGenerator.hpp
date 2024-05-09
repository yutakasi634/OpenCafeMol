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

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override { return "PositionRestraint"; }

  private:
    std::vector<std::size_t>           indices_;
    std::vector<std::array<double, 3>> positions_;
    std::vector<double>                ks_;
    std::vector<double>                v0s_;
    std::string                        ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_POSITION_RESTRAINT_FORCE_FIELD_GENERATOR_HPP
