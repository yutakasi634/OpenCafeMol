#include "ThreeSPN2CrossStackingForceFieldGenerator.hpp"

#include <algorithm>

std::unique_ptr<OpenMM::Force> ThreeSPN2CrossStackingForceFieldGenerator::generate() const noexcept
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
        "dr       = distance(d1, a3) - {id}_r0;"
        "dt3      = {id}_K_BP*(t3  - {id}_t03);"
        "dtCS     = {id}_K_CS*(tCS - {id}_t0CS);"
        "tCS      = angle(d2, d1, a3);"
        "t3       = acos(cost3lim);"
        "cost3lim = min(max(cost3, -0.99), 0.99);"
        "cost3    = sin(t1)*sin(t2)*cos(phi) - cos(t1)*cos(t2);"
        "t1       = angle(d2, d1, a1);"
        "t2       = angle(d1, a1, a2);"
        "phi      = dihedral(d2, d1, a1, a2);"
        "pi       = 3.1415926535897932385;",
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

    chbond_ff->setCutoffDistance(this->cutoff_);

    chbond_ff->addPerAcceptorParameter(fmt::format("{}_epsilon",  ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_r0",       ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_t03",      ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_t0CS",     ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_K_BP",     ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_K_CS",     ffgen_id_));
    chbond_ff->addPerAcceptorParameter(fmt::format("{}_alpha_CS", ffgen_id_));

    for (size_t i=0; i < donor_indices_vec_.size(); ++i)
    {
        const size_t d1_base_0 = donor_indices_vec_.at(i).at(0);
        const size_t d2_sugar  = donor_indices_vec_.at(i).at(1);
        const size_t d3_base_n = donor_indices_vec_.at(i).at(2);
        chbond_ff->addDonor(d1_base_0, d2_sugar, d3_base_n);
    }

    const std::string b0bp = base_pair_.first + base_pair_.second;
    for (size_t i=0; i < acceptor_indices_vec_.size(); ++i)
    {
        const size_t a1_base_p = acceptor_indices_vec_.at(i).at(0);
        const size_t a2_sugar  = acceptor_indices_vec_.at(i).at(1);
        const size_t a3_base_c = acceptor_indices_vec_.at(i).at(2);

        const auto& bpbc = base_kind_acceptors_vec_[i];

        const std::vector<double> parameters = {
            this->epsilon_CS_.at(bpbc),
            this->r0_CS_     .at(bpbc),
            this->theta3_0_  .at(b0bp),
            this->theta_CS_0_.at(bpbc),
            this->K_BP_,
            this->K_CS_,
            this->alpha_CS_
        };
        chbond_ff->addAcceptor(a1_base_p, a2_sugar, a3_base_c, parameters);
    }

    // Exclusion
    for (size_t i=0; i < donor_indices_vec_.size(); ++i)
    {
        for (size_t j=0; j < acceptor_indices_vec_.size(); ++j)
        {
            const size_t d1 = donor_indices_vec_.at(i).at(0);
            const size_t a1 = acceptor_indices_vec_.at(j).at(0);

            const auto pair = d1 < a1 ? std::make_pair(d1, a1) : std::make_pair(a1, d1);
            const auto itr  = std::find(ignore_list_.begin(), ignore_list_.end(), pair);

            if (itr != ignore_list_.end())
            {
                chbond_ff->addExclusion(i, j);
            }
        }
    }

    return chbond_ff;
}

