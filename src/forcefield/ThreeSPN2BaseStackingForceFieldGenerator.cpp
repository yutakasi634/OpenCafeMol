#include "ThreeSPN2BaseStackingForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> ThreeSPN2BaseStackingForceFieldGenerator::generate() const
{
    // Hinckley et al., J. Chem. Phys. (2013)
    std::string potential_formula = fmt::format(
        "rep + f2*attr;"                       // Eq. (8)
        " rep   = {id}_epsilon*(1 - exp(-{id}_alpha*(dr)))^2 * step(-dr);"   // Eq. (6)
        " attr  = {id}_epsilon *"
        "         (1 - exp(-{id}_alpha*(dr)))^2 * step( dr) - {id}_epsilon;" // Eq. (7)
        " dr    = distance(p2, p3) - {id}_r0;"
        " f2    = max(f*rect2, rect1);"               // Eq. (4)
        " rect1 = step(dt + pi/2) * step(pi/2 - dt);" // 1 when -pi/2 < dt < pi/2, else 0
        " rect2 = step(dt + pi)   * step(pi   - dt);" // 1 when -pi   < dt < pi,   else 0
        " f     = 1 - cos(dt)^2;"
        " dt    = {id}_K_BS * (angle(p1, p2, p3) - {id}_t0);"
        " pi    = 3.1415926535897932385;",
        fmt::arg("id", ffgen_id_));

    auto ccbond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(3, potential_formula);

    ccbond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    ccbond_ff->addPerBondParameter(fmt::format("{}_epsilon", ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_r0",      ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_t0",      ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_alpha",   ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_K_BS",    ffgen_id_));

    for (std::size_t idx=0; idx < indices_vec_.size(); ++idx)
    {
        const indices_type& triplet = indices_vec_[idx];

        std::vector<int>    particles (3);
        std::vector<double> parameters(5);

        particles [0] = triplet[0];
        particles [1] = triplet[1];
        particles [2] = triplet[2];
        parameters[0] = eps_vec_[idx];
        parameters[1] = r0s_    [idx];
        parameters[2] = theta0s_[idx];
        parameters[3] = alpha_;
        parameters[4] = K_BS_;

        ccbond_ff->addBond(particles, parameters);
    }

    return ccbond_ff;
}

