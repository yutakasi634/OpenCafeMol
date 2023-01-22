#ifndef OPEN_AICG2_PLUS_HARMONIC_ANGLE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_HARMONIC_ANGLE_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"

// [Caution!] This potential formulaw is "1/2*k*(theta - theta0)^2",
// so interaction coefficient is devided by 2. So if you use the coefficient
// base on "k*(theta  - theta0)^2" potential, you should double that value.
class HarmonicAngleForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::array<std::size_t, 3>;

  public:
    HarmonicAngleForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& v0s, const std::vector<double>& ks, const bool use_periodic)
        : indices_vec_(indices_vec), v0s_(v0s), ks_(ks), use_periodic_(use_periodic)
    {
        if(!(indices_vec.size() == v0s.size() && v0s.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] HarmonicAngleForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "v0 ("          << v0s.size()         << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                << "The number os these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        auto angle_ff = std::make_unique<OpenMM::HarmonicAngleForce>();
        angle_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& idx_triplet = indices_vec_[idx];
            angle_ff->addAngle(idx_triplet[0], idx_triplet[1], idx_triplet[2],
                               v0s_[idx], ks_[idx]);
        }

        return angle_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       v0s_;
    std::vector<double>       ks_;
    const bool                use_periodic_;
};


#endif // OPEN_AICG2_PLUS_HARMONIC_ANGLE_FORCE_FIELD_GENERATOR_HPP
