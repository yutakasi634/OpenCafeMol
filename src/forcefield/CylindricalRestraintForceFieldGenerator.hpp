#ifndef OPEN_CAFEMOL_CYLINDER_RESTRAINT_FORCEFIELD_GENERATOR_HPP
#define OPEN_CAFEMOL_CYLINDER_RESTRAINT_FORCEFIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

class CylindricalRestraintForceFieldGenerator : public ForceFieldGeneratorBase
{
  public:
    CylindricalRestraintForceFieldGenerator(
        const std::vector<std::size_t>& indices,
        const std::vector<std::array<double, 3>>& axes,
        const std::vector<std::array<double, 3>>& shifts,
        const std::vector<double>& ks, const std::vector<double>& v0s,
        const bool use_periodic)
        : indices_(indices), axes_(axes), shifts_(shifts), ks_(ks), v0s_(v0s),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("CR{}", ffid.gen()))
    {
        assert(this->indices_.size() == this->axes_.size());
        assert(this->indices_.size() == this->shifts_.size());
        assert(this->indices_.size() == this->ks_.size());
        assert(this->indices_.size() == this->v0s_.size());
        for (auto& axis : this->axes_)
        {
            const double axis_norm = std::sqrt(
                axis[0] * axis[0] +
                axis[1] * axis[1] +
                axis[2] * axis[2]);
            assert(axis_norm > 0);
            axis[0] /= axis_norm;
            axis[1] /= axis_norm;
            axis[2] /= axis_norm;
        }
        for (std::size_t idx = 0; idx < this->indices_.size(); ++idx)
        {
            const double r0ra = this->axes_[idx][0] * this->shifts_[idx][0] +
                                this->axes_[idx][1] * this->shifts_[idx][1] +
                                this->axes_[idx][2] * this->shifts_[idx][2];
            this->shifts_[idx][0] -= r0ra * this->axes_[idx][0];
            this->shifts_[idx][1] -= r0ra * this->axes_[idx][1];
            this->shifts_[idx][2] -= r0ra * this->axes_[idx][2];
        }
    }
    std::unique_ptr<OpenMM::Force> generate() const override;
    std::string name() const override { return "CylinderRestraint"; }
  private:
    std::vector<std::size_t>           indices_;
    std::vector<std::array<double, 3>> axes_;
    std::vector<std::array<double, 3>> shifts_;
    std::vector<double>                ks_;
    std::vector<double>                v0s_;
    bool                               use_periodic_;
    std::string                        ffgen_id_;
};

#endif // OPEN_CAFEMOL_CYLINDER_RESTRAINT_FORCEFIELD_GENERATOR_HPP
