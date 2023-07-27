#ifndef OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP

#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"
#include "src/util/Constants.hpp"

template<typename PotentialParameterType>
class ThreeSPN2BasePairForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 2>;
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    ThreeSPN2BasePairForceFieldGenerator(
        const std::vector<indices_type>& donor_indices_vec,
        const std::vector<indices_type>& acceptor_indices_vec,
        const std::pair<std::string, std::string>& base_pair,
        const index_pairs_type& ignore_list,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : donor_indices_vec_(donor_indices_vec),
          acceptor_indices_vec_(acceptor_indices_vec),
          base_pair_(base_pair), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
    {}

    std::unique_ptr<OpenMM::Force> generate() const override
    {
        // Hinckley et al., J. Chem. Phys. (2013)
        std::string potential_formula =
            "rep + 1/2*(1 + cos(dphi))*fdt1*fdt2*attr;" // Eq. (9)
            "rep     = epsilon*(1 - exp(-alpha_BP*dr))^2 * (1 - step(dr));"  // (1-step(dr)) = step(-dr)
            "attr    = epsilon*(1 - exp(-alpha_BP*dr))^2 * step(dr) - epsilon;"
            "fdt1    = max(f1*rect0t1, rect1t1);"
            "fdt2    = max(f2*rect0t2, rect1t2);"
            "rect1t1 = step(pi/2 + dt1) * step(pi/2 - dt1);"
            "rect1t2 = step(pi/2 + dt2) * step(pi/2 - dt2);"
            "rect0t1 = step(pi   + dt1) * step(pi   - dt1);"
            "rect0t2 = step(pi   + dt2) * step(pi   - dt2);"
            "f1      = 1 - cos(dt1)^2;"
            "f2      = 1 - cos(dt2)^2;"
            "dphi    = dihedral(d2, d1, a1, a2) - phi0;"
            "dr      = distance(d1, a1) - r0;"
            "dt1     = K_BP * (angle(d2, d1, a1) - t01);"
            "dt2     = K_BP * (angle(a2, a1, d1) - t02);"
            "pi      = 3.1415926535897932385;";

        const std::map<std::string, std::string> ff_params =
        {
            {"epsilon",  "TSPN2BP" + ffgen_id_str_ + "_epsilon"},
            {"r0",       "TSPN2BP" + ffgen_id_str_ + "_r0"},
            {"t01",      "TSPN2BP" + ffgen_id_str_ + "_t01"},
            {"t02",      "TSPN2BP" + ffgen_id_str_ + "_t02"},
            {"phi0",     "TSPN2BP" + ffgen_id_str_ + "_phi0"},
            {"alpha_BP", "TSPN2BP" + ffgen_id_str_ + "_alpha_BP"},
            {"K_BP",     "TSPN2BP" + ffgen_id_str_ + "_K_BP"},
        };

        for (const auto& param: ff_params)
        {
            potential_formula = std::regex_replace(
                potential_formula, std::regex(param.first), param.second);
        }

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
        chbond_ff->setCutoffDistance   (PotentialParameterType::cutoff);

        chbond_ff->addPerDonorParameter(ff_params.at("epsilon"));
        chbond_ff->addPerDonorParameter(ff_params.at("r0"));
        chbond_ff->addPerDonorParameter(ff_params.at("t01"));
        chbond_ff->addPerDonorParameter(ff_params.at("t02"));
        chbond_ff->addPerDonorParameter(ff_params.at("phi0"));
        chbond_ff->addPerDonorParameter(ff_params.at("alpha_BP"));
        chbond_ff->addPerDonorParameter(ff_params.at("K_BP"));

        const std::string bp_kind = base_pair_.first + base_pair_.second;
        const std::vector<double> parameters = {
            PotentialParameterType::epsilon_BP.at(bp_kind),
            PotentialParameterType::r0        .at(bp_kind),
            PotentialParameterType::theta0_1  .at(bp_kind),
            PotentialParameterType::theta0_2  .at(bp_kind),
            PotentialParameterType::phi0      .at(bp_kind),
            PotentialParameterType::alpha_BP,
            PotentialParameterType::K_BP,
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

    std::string name() const noexcept
    {
        return PotentialParameterType::name+"BasePair "
               "(" + base_pair_.first + "-" + base_pair_.second + ")";
    }

  private:
    std::vector<indices_type> donor_indices_vec_;
    std::vector<indices_type> acceptor_indices_vec_;
    std::pair<std::string, std::string> base_pair_;
    index_pairs_type          ignore_list_;
    bool                      use_periodic_;
    std::string               ffgen_id_str_;
};


// ----------------------------------------------------------------------------
// 3SPN2 parameters
//
// The parameters are derived from the following literature
// - W. Lu, C. Bueno, N. P. Schafer, J. Moller, S. Jin, X. Chen, M. Chen, X. Gu,
//   A. Davtyan, J. J. de Pablo, and P. G. Wolynes, PLOS Comput. Biol. (2021)
//
//   - Table 10
//
// Also see
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)

struct ThreeSPN2BasePairPotentialParameter
{
    inline static double cutoff   = 18.0 * OpenMM::NmPerAngstrom;  // [nm]

    inline static double alpha_BP =  2.0 / OpenMM::NmPerAngstrom;  // [1/nm]
    inline static double K_BP     = 12.0;

    inline static const std::map<std::string, double> epsilon_BP = { // [kJ/mol]
        {"AT", 16.73},
        {"TA", 16.73},
        {"GC", 21.18},
        {"CG", 21.18}
    };

    inline static const std::map<std::string, double> r0 = { // [nm]
        {"AT", 5.941 * OpenMM::NmPerAngstrom},
        {"TA", 5.941 * OpenMM::NmPerAngstrom},
        {"GC", 5.530 * OpenMM::NmPerAngstrom},
        {"GC", 5.530 * OpenMM::NmPerAngstrom}
    };

    inline static const std::map<std::string, double> theta0_1 = { // [radian]
        {"AT", 156.54 * OpenMM::RadiansPerDegree},
        {"TA", 135.78 * OpenMM::RadiansPerDegree},
        {"GC", 159.81 * OpenMM::RadiansPerDegree},
        {"GC", 141.16 * OpenMM::RadiansPerDegree}
    };

    inline static const std::map<std::string, double> theta0_2 = { // [radian]
        {"AT", 135.78 * OpenMM::RadiansPerDegree},
        {"TA", 156.54 * OpenMM::RadiansPerDegree},
        {"GC", 141.16 * OpenMM::RadiansPerDegree},
        {"GC", 159.81 * OpenMM::RadiansPerDegree}
    };

    inline static const std::map<std::string, double> phi0 = { // [radian]
        {"AT", -38.35 * OpenMM::RadiansPerDegree},
        {"TA", -38.35 * OpenMM::RadiansPerDegree},
        {"GC", -42.98 * OpenMM::RadiansPerDegree},
        {"GC", -42.98 * OpenMM::RadiansPerDegree}
    };

    inline static const std::string name = "3SPN2";
};

// ----------------------------------------------------------------------------
// 3SPN2 parameters
//
// The parameters are derived from the following literature
// - W. Lu, C. Bueno, N. P. Schafer, J. Moller, S. Jin, X. Chen, M. Chen, X. Gu,
//   A. Davtyan, J. J. de Pablo, and P. G. Wolynes, PLOS Comput. Biol. (2021)
//
//   - Table 10
//
// Also see
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// - G. S. Freeman, D. M. Hinckley, J. P. Lequieu, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2014)

struct ThreeSPN2CBasePairPotentialParameter
{
    inline static const double cutoff       = 18.0 * OpenMM::NmPerAngstrom;  // [nm]

    inline static const double alpha_BP     =  2.0 / OpenMM::NmPerAngstrom;  // [1/nm]
    inline static const double K_BP         = 12.0;

    inline static const std::map<std::string, double> epsilon_BP = { // [kJ/mol]
        {"AT", 14.41},
        {"TA", 14.41},
        {"GC", 18.24},
        {"CG", 18.24}
    };

    inline static const std::map<std::string, double> r0 = { // [nm]
        {"AT", 5.82 * OpenMM::NmPerAngstrom},
        {"TA", 5.82 * OpenMM::NmPerAngstrom},
        {"GC", 5.52 * OpenMM::NmPerAngstrom},
        {"GC", 5.52 * OpenMM::NmPerAngstrom}
    };

    inline static const std::map<std::string, double> theta0_1 = { // [radian]
        {"AT", 153.17 * OpenMM::RadiansPerDegree},
        {"TA", 133.51 * OpenMM::RadiansPerDegree},
        {"GC", 159.50 * OpenMM::RadiansPerDegree},
        {"GC", 138.08 * OpenMM::RadiansPerDegree}
    };

    inline static const std::map<std::string, double> theta0_2 = { // [radian]
        {"AT", 133.51 * OpenMM::RadiansPerDegree},
        {"TA", 153.17 * OpenMM::RadiansPerDegree},
        {"GC", 138.08 * OpenMM::RadiansPerDegree},
        {"GC", 159.50 * OpenMM::RadiansPerDegree}
    };

    inline static const std::map<std::string, double> phi0 = { // [radian]
        {"AT", -38.18 * OpenMM::RadiansPerDegree},
        {"TA", -38.18 * OpenMM::RadiansPerDegree},
        {"GC", -35.75 * OpenMM::RadiansPerDegree},
        {"GC", -35.75 * OpenMM::RadiansPerDegree}
    };
    inline static const std::string name = "3SPN2C";
};

#endif // OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP
