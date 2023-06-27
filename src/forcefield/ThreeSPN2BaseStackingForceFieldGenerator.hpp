#ifndef OPEN_AICG2_PLUS_3SPN2_BASE_STACKING_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BASE_STACKING_FORCE_FIELD_GENERATOR_HPP

#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>

#include "src/util/Constants.hpp"

class ThreeSPN2BaseStackingForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::array<std::size_t, 3>;

  public:
    ThreeSPN2BaseStackingForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& eps_vec, const std::vector<double>& r0s,
        const std::vector<double>& theta0, const double& alpha,
        const double& K_BS,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : indices_vec_(indices_vec), eps_vec_(eps_vec), r0s_(r0s),
          theta0s_(theta0), alpha_(alpha), K_BS_(K_BS),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
    {
        if(!(indices_vec.size() == eps_vec.size() && eps_vec.size() == r0s.size() &&
             r0s.size() == theta0.size()))
        {
            std::ostringstream oss;
            oss << "[error] ThreeSPN2BaseStackingForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "eps_vec ("     << eps_vec.size()     << "), "
                   "r0s ("         << r0s.size()         << "),"
                   "theta0 ("      << theta0.size()      << ") is not matched "
               << "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    };

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        // Hinckley et al., J. Chem. Phys. (2013)
        std::string potential_formula = "energy;"
            " energy = rep + f2*attr;"                                          // Eq. (8)
            " rep    = epsilon*(1 - exp(-alpha*(dr)))^2 * step(-dr);"           // Eq. (6)
            " attr   = epsilon*(1 - exp(-alpha*(dr)))^2 * step( dr) - epsilon;" // Eq. (7)
            " dr     = distance(p2, p3) - r0;"
            " f2     = max(f*rect2, rect1);"               // Eq. (4)
            " rect1  = step(dt + pi/2) * step(pi/2 - dt);" // 1 when -pi/2 < dt < pi/2, else 0
            " rect2  = step(dt + pi)   * step(pi   - dt);" // 1 when -pi   < dt < pi,   else 0
            " f      = 1 - cos(dt)^2;"
            " dt     = K_BS * (angle(p1, p2, p3) - t0);"
            " pi     = 3.1415926535897932385;";

        const std::map<std::string, std::string> ff_params =
        {
            {"epsilon", "TSPN2BS" + ffgen_id_str_ + "_epsilon"},
            {"r0",      "TSPN2BS" + ffgen_id_str_ + "_r0"},
            {"t0",      "TSPN2BS" + ffgen_id_str_ + "_t0"},
            {"alpha",   "TSPN2BS" + ffgen_id_str_ + "_alpha"},
            {"K_BS",    "TSPN2BS" + ffgen_id_str_ + "_K_BS"},
        };

        for (auto itr=ff_params.begin(); itr!=ff_params.end(); ++itr)
        {
            potential_formula = std::regex_replace(
                potential_formula, std::regex(itr->first), itr->second);
        }
        auto ccbond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(3, potential_formula);

        ccbond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        ccbond_ff->addPerBondParameter(ff_params.at("epsilon"));
        ccbond_ff->addPerBondParameter(ff_params.at("r0"));
        ccbond_ff->addPerBondParameter(ff_params.at("t0"));
        ccbond_ff->addPerBondParameter(ff_params.at("alpha"));
        ccbond_ff->addPerBondParameter(ff_params.at("K_BS"));

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

    std::string name() const noexcept {return "3SPN2BaseStacking";};

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       eps_vec_;
    std::vector<double>       r0s_;
    std::vector<double>       theta0s_;
    const double              alpha_;
    const double              K_BS_;
    const bool                use_periodic_;
    const std::string         ffgen_id_str_;
};

// ----------------------------------------------------------------------------
// 3SPN2 parameters
//
// The parameters are derived from the following literature
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
//   - TABLE IV:  K_BS, alph_BS
//   - TABLE VII: epsilon_BS, r0_BS, theta0_BS

struct ThreeSPN2BaseStackingPotentialParameter
{
    double K_BS         = 6.0;
    double alpha_BS     = 3.0; // [1/angstrom]

    std::map<std::string, double> epsilon_BS = { // [kJ/mol]
        {"AA", 14.39}, {"AT", 14.34}, {"AG", 13.25}, {"AC", 14.51},
        {"TA", 10.37}, {"TT", 13.36}, {"TG", 10.34}, {"TC", 12.89},
        {"GA", 14.81}, {"GT", 15.57}, {"GG", 14.93}, {"GC", 15.39},
        {"CA", 11.42}, {"CT", 12.79}, {"CG", 10.52}, {"CC", 13.24},
    };

    std::map<std::string, double> r0_BS = { // [angstrom]
        {"AA", 3.716}, {"AT", 3.675}, {"AG", 3.827}, {"AC", 3.744},
        {"TA", 4.238}, {"TT", 3.984}, {"TG", 4.416}, {"TC", 4.141},
        {"GA", 3.576}, {"GT", 3.598}, {"GG", 3.664}, {"GC", 3.635},
        {"CA", 4.038}, {"CT", 3.798}, {"CG", 4.208}, {"CC", 3.935},
    };

    std::map<std::string, double> theta0_BS = { // [degree]
        {"AA", 101.15}, {"AT", 85.94}, {"AG", 105.26}, {"AC", 89.00},
        {"TA", 101.59}, {"TT", 89.50}, {"TG", 104.31}, {"TC", 91.28},
        {"GA", 100.89}, {"GT", 84.83}, {"GG", 105.48}, {"GC", 88.28},
        {"CA", 106.49}, {"CT", 93.31}, {"CG", 109.54}, {"CC", 95.46},
    };
};

// ----------------------------------------------------------------------------
// 3SPN2.C parameters
//
// The parameters are derived from the following literature
//
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
//   - TABLE IV: K_BS, alpha_BS
//
// - G. S. Freeman, D. M. Hinckley, J. P. Lequieu, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2014)
//
//   - TABLE VII: epsilon_BS
//
// - W. Lu, C. Bueno, N. P. Schafer, J. Moller, S. Jin, X. Chen, M. Chen, X. Gu,
//   A. Davtyan, J. J. de Pablo, and P. G. Wolynes, PLOS Comput. Biol. (2021)
//
//   - Table 9: r0_BS, theta0_BS

struct ThreeSPN2CBaseStackingPotentialParameter
{
    double K_BS         = 6.0;
    double alpha_BS     = 3.0; // [1/angstrom]

    std::map<std::string, double> epsilon_BS = { // [kJ/mol]
        {"AA", 13.82},  {"AT", 15.05}, {"AG",  13.32}, {"AC",  15.82},
        {"TA",  9.15},  {"TT", 12.44}, {"TG",   9.58}, {"TC",  13.11},
        {"GA", 13.76},  {"GT", 14.59}, {"GG",  14.77}, {"GC",  15.17},
        {"CA",  9.25},  {"CT", 12.42}, {"CG",   8.83}, {"CC",  14.01},
    };

    std::map<std::string, double> r0_BS = { // [angstrom]
        {"AA", 3.58},   {"AT", 3.56},  {"AG",   3.85}, {"AC",   3.45},
        {"TA", 4.15},   {"TT", 3.93},  {"TG",   4.32}, {"TC",   3.87},
        {"GA", 3.51},   {"GT", 3.47},  {"GG",   3.67}, {"GC",   3.42},
        {"CA", 4.15},   {"CT", 3.99},  {"CG",   4.34}, {"CC",   3.84},
    };

    std::map<std::string, double> theta0_BS = { // [degree]
        {"AA", 100.13}, {"AT", 90.48}, {"AG", 104.39}, {"AC",  93.23},
        {"TA", 102.59}, {"TT", 93.32}, {"TG", 103.70}, {"TC",  94.55},
        {"GA",  95.45}, {"GT", 87.63}, {"GG", 106.36}, {"GC",  83.12},
        {"CA", 102.69}, {"CT", 96.05}, {"CG", 100.46}, {"CC", 100.68},
    };
};

#endif // OPEN_AICG2_PLUS_3SPN2_BASE_STACKING_FORCE_FIELD_GENERATOR_HPP
