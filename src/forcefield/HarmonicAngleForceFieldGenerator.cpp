#include "HarmonicAngleForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> HarmonicAngleForceFieldGenerator::generate() const
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
