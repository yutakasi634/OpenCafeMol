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

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept override { return "Pulling"; }

  private:
    std::vector<parameter_type> params_;
    bool                        use_periodic_;
    std::string                 ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_PULLING_FORCE_FIELD_GENERATOR_HPP
