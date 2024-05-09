#include "ThreeSPN2BasePairForceFieldGenerator.hpp"

#include <algorithm>

std::unique_ptr<OpenMM::Force> ThreeSPN2BasePairForceFieldGenerator::generate() const
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
        "dphi    = dihedral(d2, d1, a1, a2) - {id}_phi0;"
        "dr      = distance(d1, a1) - {id}_r0;"
        "dt1     = {id}_K_BP * (angle(d2, d1, a1) - {id}_t01);"
        "dt2     = {id}_K_BP * (angle(a2, a1, d1) - {id}_t02);"
        "pi      = 3.1415926535897932385;",
        fmt::arg("id", ffgen_id_));

    auto chbond_ff = std::make_unique<OpenMM::CustomHbondForce>(potential_formula);

    // set cutoff
    if(use_periodic_)
    {
        chbond_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
    }
    else
    {
        chbond_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffNonPeriodic);
    }

    // set cutoff
    chbond_ff->setCutoffDistance   (this->cutoff_);

    chbond_ff->addPerDonorParameter(fmt::format("{}_epsilon",  ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_r0",       ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_t01",      ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_t02",      ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_phi0",     ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_alpha_BP", ffgen_id_));
    chbond_ff->addPerDonorParameter(fmt::format("{}_K_BP",     ffgen_id_));

    const std::string bp_kind = base_pair_.first + base_pair_.second;
    const std::vector<double> parameters = {
        this->epsilon_BP_.at(bp_kind),
        this->r0_        .at(bp_kind),
        this->theta0_1_  .at(bp_kind),
        this->theta0_2_  .at(bp_kind),
        this->phi0_      .at(bp_kind),
        this->alpha_BP_,
        this->K_BP_
    };

    for (const auto& donor_particles: donor_indices_vec_)
    {
        const size_t d1_base  = donor_particles.at(0);
        const size_t d2_sugar = donor_particles.at(1);
        chbond_ff->addDonor(d1_base, d2_sugar, -1, parameters);
    }

    for (const auto& acceptor_particles: acceptor_indices_vec_)
    {
        const size_t a1_base  = acceptor_particles.at(0);
        const size_t a2_sugar = acceptor_particles.at(1);
        chbond_ff->addAcceptor(a1_base, a2_sugar, -1);
    }

    // set exclusion list
    for (const auto& donor_particles: donor_indices_vec_)
    {
        const size_t d_idx = &donor_particles - &donor_indices_vec_[0];
        const size_t d1    = donor_particles.at(0); // base particle

        for (const auto& acceptor_particles: acceptor_indices_vec_)
        {
            const size_t a_idx = &acceptor_particles - &acceptor_indices_vec_[0];
            const size_t a1    = acceptor_particles.at(0); // base particle

            // Find whether base-particles pair in this donor-acceptor is included in ignore_list
            const auto pair = d1 < a1 ? std::make_pair(d1, a1): std::make_pair(a1, d1);
            const auto itr  = std::find(ignore_list_.begin(), ignore_list_.end(), pair);

            if (itr != ignore_list_.end())
            {
                chbond_ff->addExclusion(d_idx, a_idx);
            }
        }
    }

    return chbond_ff;
}

