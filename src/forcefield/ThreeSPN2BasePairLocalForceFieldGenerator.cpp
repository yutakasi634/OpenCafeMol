#include "ThreeSPN2BasePairLocalForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> ThreeSPN2BasePairLocalForceFieldGenerator::generate() const
{
    // Hinckley et al., J. Chem. Phys. (2013)
    std::string potential_formula = fmt::format(
        "rep + 1/2*(1 + cos(dphi))*fdt1*fdt2*attr;" // Eq. (9)
        "rep     = {id}_epsilon*(1 - exp(-{id}_alpha_BP*dr))^2 * (1 - step(dr));"  // (1-step(dr)) = step(-dr)
        "attr    = {id}_epsilon *"
        "          (1 - exp(-{id}_alpha_BP*dr))^2 * step(dr) - {id}_epsilon;"
        "fdt1    = max(f1*rect0t1, rect1t1);"
        "fdt2    = max(f2*rect0t2, rect1t2);"
        "rect1t1 = step(pi/2 + dt1) * step(pi/2 - dt1);"
        "rect1t2 = step(pi/2 + dt2) * step(pi/2 - dt2);"
        "rect0t1 = step(pi   + dt1) * step(pi   - dt1);"
        "rect0t2 = step(pi   + dt2) * step(pi   - dt2);"
        "f1      = 1 - cos(dt1)^2;"
        "f2      = 1 - cos(dt2)^2;"
        "dphi    = dihedral(p1, p2, p3, p4) - {id}_phi0;"
        "dr      = distance(p2, p3) - {id}_r0;"
        "dt1     = {id}_K_BP * (angle(p1, p2, p3) - {id}_t01);"
        "dt2     = {id}_K_BP * (angle(p4, p3, p2) - {id}_t02);"
        "pi      = 3.1415926535897932385;",
        fmt::arg("id", ffgen_id_));

    auto ccbond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(4, potential_formula);

    ccbond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    ccbond_ff->addPerBondParameter(fmt::format("{}_epsilon",  ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_r0",       ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_t01",      ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_t02",      ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_phi0",     ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_alpha_BP", ffgen_id_));
    ccbond_ff->addPerBondParameter(fmt::format("{}_K_BP",     ffgen_id_));

    for(const auto& indices : indices_vec_)
    {
        const std::vector<int> particles = {
            static_cast<int>(indices[0]),
            static_cast<int>(indices[1]),
            static_cast<int>(indices[2]),
            static_cast<int>(indices[3]),
        };

        const std::string bp_kind = base_pair_.first + base_pair_.second;
        const std::vector<double> parameters = {
            this->epsilon_BP_.at(bp_kind),
            this->r0_        .at(bp_kind),
            this->theta0_1_  .at(bp_kind),
            this->theta0_2_  .at(bp_kind),
            this->phi0_      .at(bp_kind),
            this->alpha_BP_,
            this->K_BP_,
        };

        ccbond_ff->addBond(particles, parameters);
    }

    return ccbond_ff;
}

