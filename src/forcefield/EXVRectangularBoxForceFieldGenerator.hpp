#ifndef OPEN_AICG2_PLUS_EXV_RECTANGULAR_BOX_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_EXV_RECTANGULAR_BOX_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class EXVRectangularBoxForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    EXVRectangularBoxForceFieldGenerator(
        const double eps,
        const std::array<double, 3>& box_lower, const std::array<double, 3>& box_upper,
        const std::vector<std::size_t>& indices, const std::vector<double>& radii,
        const bool use_periodic);

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override { return "EXVRectangularBox"; }

  private:
    double                   eps_;
    std::array<double, 3>    box_lower_;
    std::array<double, 3>    box_upper_;
    std::vector<std::size_t> indices_;
    std::vector<double>      radii_;
    bool                     use_periodic_;
    std::string              ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_EXV_RECTANGULAR_BOX_FORCE_FIELD_GENERATOR_HPP
