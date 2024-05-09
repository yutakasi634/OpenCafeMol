#include "ThreeSPN2CrossStackingLocalForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> ThreeSPN2CrossStackingLocalForceFieldGenerator::generate() const
{
    // Hinckley et al., J. Chem. Phys. (2013)
    std::string potential_formula = fmt::format(
        "fdt3*fdtCS*attr/2;"  // Eq. (10)
        "attr     = {id}_epsilon*"
        "           (1 - exp(-{id}_alpha_CS*dr))^2 * step(dr) - {id}_epsilon;"
        "fdt3     = max(f1*rect0t3,  rect1t3);"
        "fdtCS    = max(f2*rect0tCS, rect1tCS);"
        "rect0t3  = step(pi   + dt3)  * step(pi   - dt3);"
        "rect0tCS = step(pi   + dtCS) * step(pi   - dtCS);"
        "rect1t3  = step(pi/2 + dt3)  * step(pi/2 - dt3);"
        "rect1tCS = step(pi/2 + dtCS) * step(pi/2 - dtCS);"
        "f1       = 1 - cos(dt3)^2;"
        "f2       = 1 - cos(dtCS)^2;"
        "dr       = distance(p2, p6) - {id}_r0;"
        "dt3      = {id}_K_BP*(t3  - {id}_t03);"
        "dtCS     = {id}_K_CS*(tCS - {id}_t0CS);"
        "tCS      = angle(p1, p2, p6);"
        "t3       = acos(cost3lim);"
        "cost3lim = min(max(cost3, -0.99), 0.99);"
        "cost3    = sin(t1)*sin(t2)*cos(phi) - cos(t1)*cos(t2);"
        "t1       = angle(p1, p2, p3);"
        "t2       = angle(p2, p3, p4);"
        "phi      = dihedral(p1, p2, p3, p4);"
        "pi       = 3.1415926535897932385;",
        fmt::arg("id", ffgen_id_));

    auto ccbond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(6, potential_formula);

    ccbond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    ccbond_ff->addPerBondParameter(fmt::format("{}_epsilon",  ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_r0",       ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_t03",      ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_t0CS",     ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_K_BP",     ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_K_CS",     ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_alpha_CS", ffgen_id_));

    for (std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const indices_type& indices = indices_vec_[idx];
        const std::vector<int> particles = {
            static_cast<int>(indices[0]),
            static_cast<int>(indices[1]),
            static_cast<int>(indices[2]),
            static_cast<int>(indices[3]),
            static_cast<int>(indices[4]),
            static_cast<int>(indices[5]),
        };

        const auto&       bpbc = base_kind_vec_[idx];
        const std::string b0bp = base_pair_.first + base_pair_.second;
        const std::vector<double> parameters = {
            this->epsilon_CS_.at(bpbc),
            this->r0_CS_     .at(bpbc),
            this->theta3_0_  .at(b0bp),
            this->theta_CS_0_.at(bpbc),
            this->K_BP_,
            this->K_CS_,
            this->alpha_CS_
        };

        ccbond_ff->addBond(particles, parameters);
    }

    return ccbond_ff;
}
