#ifndef OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP

#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "src/util/Constants.hpp"

class ProteinDNANonSpecificForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_dna_type     = std::array<std::size_t, 2>;
    using indices_protein_type = std::array<std::size_t, 3>;

  public:
    ProteinDNANonSpecificForceFieldGenerator(
        const std::vector<indices_dna_type>& indices_dna,
        const std::vector<indices_protein_type>& indices_protein,
        const double sigma, const double delta,
        const std::vector<double> ks, const std::vector<double> r0s,
        const std::vector<double> theta0s, const std::vector<double> phi0s,
        const bool use_periodic)
        : indices_dna_(indices_dna), indices_protein_(indices_protein),
          sigma_(sigma), delta_(delta), ks_(ks), r0s_(r0s),
          theta0s_(theta0s), phi0s_(phi0s),
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
            "- {id}_k * kT * kcal2kJ * f_r * g_theta * g_phi;"
            "f_r       = exp(-dr^2/ (2 * {id}_sigma * {id}_sigma));"
            "g_theta   = max(g1*rect0t, rect1t);"
            "g_phi     = max(g2*rect0p, rect1p);"
            "rect0t    = step(pi   + dtheta) * step(pi   - dtheta);"  // 1 when -pi   < dt < pi  , else 0
            "rect1t    = step(pi/2 + dtheta) * step(pi/2 - dtheta);"  // 1 when -pi/2 < dt < pi/2, else 0
            "rect0p    = step(pi   + dphi)   * step(pi   - dphi);"    // 1 when -pi   < dt < pi  , else 0
            "rect1p    = step(pi/2 + dphi)   * step(pi/2 - dphi);"    // 1 when -pi/2 < dt < pi/2, else 0
            "g1        = 1 - cos(dtheta)^2;"
            "g2        = 1 - cos(dphi)^2;"
            "dr        = distance(p2, p4) - {id}_r0;"
            "dtheta    = K * (theta - {id}_theta0);"
            "dphi      = K * (phi   - {id}_phi0);"
            "K         = pi/(2 * {id}_delta);"
            "phi       = angle(p2, p4, p5);"
            "theta     = acos(cost1lim);"           // angle between vecots CA_N->CA_C and CA->P
            "cost1lim  = min(max(cost1, -1), 1);"
            "cost1     = sin(t1)*sin(t2)*cos(phi1) - cos(t1)*cos(t2);"
            "t1        = angle(p3, p1, p4);"
            "t2        = angle(p1, p4, p2);"
            "phi1      = dihedral(p3, p1, p4, p2);" // p3->p1 (vector CA_N -> CA_C), p2->p4 (vector CA->P)
            "kcal2kJ   = 4.184;"
            "kT        = 0.593;"
            "pi        = 3.1415926535897932385;",
            fmt::arg("id", ffgen_id_));

        auto ccbond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(5, potential_formula);

        ccbond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        ccbond_ff->addPerBondParameter(fmt::format("{}_sigma",  ffgen_id_));
        ccbond_ff->addPerBondParameter(fmt::format("{}_delta",  ffgen_id_));
        ccbond_ff->addPerBondParameter(fmt::format("{}_k",      ffgen_id_));
        ccbond_ff->addPerBondParameter(fmt::format("{}_r0",     ffgen_id_));
        ccbond_ff->addPerBondParameter(fmt::format("{}_theta0", ffgen_id_));
        ccbond_ff->addPerBondParameter(fmt::format("{}_phi0",   ffgen_id_));

        for (const auto& idxs_pro: indices_protein_)
        {
            const size_t ipro = &idxs_pro - &indices_protein_[0];
            
            for(const auto& idxs_dna: indices_dna_)
            {
                const std::vector<int> particles = {
                    static_cast<int>(idxs_pro[2]), // CA_C   (Protein)
                    static_cast<int>(idxs_pro[0]), // CA     (Protein)
                    static_cast<int>(idxs_pro[1]), // CA_N   (Protein)
                    static_cast<int>(idxs_dna[0]), // Phos   (DNA)
                    static_cast<int>(idxs_dna[1]), // Sugar3 (DNA)
                };

                std::cout << "KJPerKcal = " << OpenMM::KJPerKcal << std::endl;
                std::cout << "KcalPerKJ = " << OpenMM::KcalPerKJ << std::endl;
                const std::vector<double> parameters = {
                    sigma_,
                    delta_,
                    ks_[ipro],
                    r0s_[ipro],
                    theta0s_[ipro],
                    phi0s_[ipro],
                };
                ccbond_ff->addBond(particles, parameters);
            }
        }

        return ccbond_ff;
    }

    std::string name() const noexcept { return "PDNS"; }

  private:
    std::vector<indices_dna_type> indices_dna_;
    std::vector<indices_protein_type> indices_protein_;
    double sigma_;
    double delta_;
    std::vector<double> ks_;
    std::vector<double> r0s_;
    std::vector<double> theta0s_;
    std::vector<double> phi0s_;
    bool                use_periodic_;
    std::string         ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
