#ifndef OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP

#include "src/util/Constants.hpp"
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>


class ProteinDNANonSpecificForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_dna_type     = std::array<std::size_t, 2>;
    using indices_protein_type = std::array<std::size_t, 3>;

  public:
    ProteinDNANonSpecificForceFieldGenerator(
        const std::vector<indices_dna_type>& indices_dna,
        const std::vector<indices_protein_type>& indices_protein,
        const double sigma, const double delta, const double cutoff_ratio,
        const std::vector<double> ks, const std::vector<double> r0s,
        const std::vector<double> theta0s, const std::vector<double> phi0s,
        const bool use_periodic)
        : indices_dna_(indices_dna), indices_protein_(indices_protein),
          sigma_(sigma), delta_(delta), cutoff_ratio_(cutoff_ratio),
          ks_(ks), r0s_(r0s), theta0s_(theta0s), phi0s_(phi0s),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("PDNS{}", ffid.gen()))
    {
        if(!(indices_protein.size() == ks.size() && ks.size() == r0s.size() &&
              r0s.size() == theta0s.size() && theta0s.size() == phi0s.size()))
        {
            std::ostringstream oss;
            oss << "[error] ProteinDNANonSpecificForceFieldGenerator: "
                   "parameter number of "
                   "indices_dna (" << indices_dna.size() << "), "
                   "k ("           << ks.size()          << "),"
                   "r0 ("          << r0s.size()         << "),"
                   "theta0 ("      << theta0s.size()     << "), and "
                   "phi0 ("        << phi0s.size()       << ") is not matched. "
                << "The number os these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override
    {
        // T. Niina et al., PloS Comp. Biol (2017)
        // U(r, theta, phi) = k*f(r)*g(theta)*g(phi)
        // 
        // where,
        //    f(r)   = exp(-(r-r0)^2 / 2signam^2)
        //    g(phi) = 1                              ...         |phi - phi0| <  delta
        //             1 - cos^2(pi(phi-phi0)/2delta) ... delta < |phi - phi0| < 2delta
        //             0                              ... otherwise

        /*
        //  PC          S5'
        //    o         o--o B
        //    A\ P   D /
        // vP | o --> o
        //    |/     `-\ 
        //    o    phi  o--o B
        //  PN          S3'
        */

        std::string potential_formula = fmt::format(
            "- {id}_k * f_r * g_theta * g_phi;"
            "f_r       = exp(-dr^2/ (2 * {id}_sigma * {id}_sigma));"
            "g_theta   = max(g1*rect0t, rect1t);"
            "g_phi     = max(g2*rect0p, rect1p);"
            "rect0t    = step(pi   + dtheta) * step(pi   - dtheta);"  // 1 when -pi   < dt < pi  , else 0
            "rect1t    = step(pi/2 + dtheta) * step(pi/2 - dtheta);"  // 1 when -pi/2 < dt < pi/2, else 0
            "rect0p    = step(pi   + dphi)   * step(pi   - dphi);"    // 1 when -pi   < dt < pi  , else 0
            "rect1p    = step(pi/2 + dphi)   * step(pi/2 - dphi);"    // 1 when -pi/2 < dt < pi/2, else 0
            "g1        = 1 - cos(dtheta)^2;"
            "g2        = 1 - cos(dphi)^2;"
            "dr        = distance(d2, a1) - {id}_r0;"
            "dtheta    = K * (theta - {id}_theta0);"
            "dphi      = K * (phi   - {id}_phi0);"
            "K         = pi/(2 * {id}_delta);"
            "phi       = angle(d2, a1, a2);"
            "theta     = acos(cost1lim);"// angle between vecots CA_N->CA_C and CA->P
            "cost1lim  = min(max(cost1, -0.99999), 0.99999);" // clliping to avoid the acos singularity due to cost1 taking values of -1 or 1.
            "cost1     = sin(t1)*sin(t2)*cos(phi1) - cos(t1)*cos(t2);"
            "t1        = angle(d1, d3, a1);"
            "t2        = angle(d3, a1, d2);"
            "phi1      = dihedral(d1, d3, a1, d2);" // d1->d3 (vector CA_N -> CA_C), d2->a1 (vector CA->P)
            "pi        = 3.1415926535897932385;",
            fmt::arg("id", ffgen_id_));

        auto chbond_ff = std::make_unique<OpenMM::CustomHbondForce>(potential_formula);

        chbond_ff->addPerDonorParameter(fmt::format("{}_sigma",  ffgen_id_));
        chbond_ff->addPerDonorParameter(fmt::format("{}_delta",  ffgen_id_));
        chbond_ff->addPerDonorParameter(fmt::format("{}_k",      ffgen_id_));
        chbond_ff->addPerDonorParameter(fmt::format("{}_r0",     ffgen_id_));
        chbond_ff->addPerDonorParameter(fmt::format("{}_theta0", ffgen_id_));
        chbond_ff->addPerDonorParameter(fmt::format("{}_phi0",   ffgen_id_));

        for (const auto& donor_particles: indices_protein_)
        {
            const size_t idx = &donor_particles - &indices_protein_[0];

            const size_t d1_calpha_n = donor_particles.at(1);
            const size_t d2_calpha   = donor_particles.at(0);
            const size_t d3_calpha_c = donor_particles.at(2);

            const std::vector<double> parameters = {
                sigma_,
                delta_,
                ks_[idx],
                r0s_[idx],
                theta0s_[idx],
                phi0s_[idx],
            };
            chbond_ff->addDonor(d1_calpha_n, d2_calpha, d3_calpha_c, parameters);
        }

        for (const auto& acceptor_particles: indices_dna_)
        {
            const size_t a1_phos   = acceptor_particles.at(0);
            const size_t a2_sugar3 = acceptor_particles.at(1);
            chbond_ff->addAcceptor(a1_phos, a2_sugar3, -1);
        }

        // set cutoff
        if(use_periodic_)
        {
            chbond_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
        }
        else
        {
            chbond_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffNonPeriodic);
        }

        double max_cutoff_length = 0.0;
        for (const auto& r0: r0s_)
        {
            max_cutoff_length = std::max(max_cutoff_length, r0 + cutoff_ratio_ * sigma_);
        }
        chbond_ff->setCutoffDistance(max_cutoff_length);

        std::cerr << "    PDNS                          : cutoff distance is "
                  << max_cutoff_length << " nm" << std::endl;

        return chbond_ff;
    }

    std::string name() const noexcept { return "PDNS"; }

  private:
    std::vector<indices_dna_type> indices_dna_;
    std::vector<indices_protein_type> indices_protein_;
    double sigma_;
    double delta_;
    double cutoff_ratio_;
    std::vector<double> ks_;
    std::vector<double> r0s_;
    std::vector<double> theta0s_;
    std::vector<double> phi0s_;
    bool                use_periodic_;
    std::string         ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
