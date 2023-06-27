#ifndef OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP

#include <cmath>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>

#include "src/util/Constants.hpp"

template<typename realT>
class ThreeSPN2CrossStackingParameterList
{
  using real_type = realT;

  public:
    template<typename ParameterSet>
    ThreeSPN2CrossStackingParameterList(ParameterSet param_set,
        const std::string& b0bp, const std::string& bpbc)
    : cutoff_    (param_set.cutoff),
      alpha_CS_  (param_set.alpha_CS),
      K_CS_      (param_set.K_CS),
      K_BP_      (param_set.K_BP),
      epsilon_CS_(param_set.epsilon_CS.at(bpbc)),
      r0_CS_     (param_set.r0_CS     .at(bpbc)),
      theta3_0_  (param_set.theta3_0  .at(b0bp)),
      theta_CS_0_(param_set.theta_CS_0.at(bpbc)),
      b0bp_  (b0bp),
      bpbc_  (bpbc)
    {}

    real_type cutoff()     const noexcept {return cutoff_;}
    real_type alpha_CS()   const noexcept {return alpha_CS_;}
    real_type K_CS()       const noexcept {return K_CS_;}
    real_type K_BP()       const noexcept {return K_BP_;}
    real_type epsilon_CS() const noexcept {return epsilon_CS_;}
    real_type r0_CS()      const noexcept {return r0_CS_;}
    real_type theta3_0()   const noexcept {return theta3_0_;}
    real_type theta_CS_0() const noexcept {return theta_CS_0_;}
    std::string getBasePairBases()      const noexcept {return b0bp_;}
    std::string getCrossStackingBases() const noexcept {return bpbc_;}

  private:
    const real_type  cutoff_;
    const real_type  alpha_CS_;
    const real_type  K_CS_;
    const real_type  K_BP_;
    const real_type  epsilon_CS_;
    const real_type  r0_CS_;
    const real_type  theta3_0_;
    const real_type  theta_CS_0_;
    const std::string b0bp_;
    const std::string bpbc_;
};


class ThreeSPN2CrossStackingForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 3>;
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;
    using parameter_list   = ThreeSPN2CrossStackingParameterList<double>;

  public:
    ThreeSPN2CrossStackingForceFieldGenerator(
        const std::vector<indices_type>& donor_indices_vec,
        const std::vector<indices_type>& acceptor_indices_vec,
        const index_pairs_type& ignore_list,
        const std::vector<parameter_list>& parameters_vec,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : donor_indices_vec_(donor_indices_vec), acceptor_indices_vec_(acceptor_indices_vec),
          ignore_list_(ignore_list), parameters_vec_(parameters_vec),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
    {
        if(!(acceptor_indices_vec.size() == parameters_vec.size()))
        {
            std::ostringstream oss;
            oss << "[error] ThreeSPN2CrossStackingForceFieldGenerator: "
                   "parameter number of "
                   "acceptor_indices_vec (" << acceptor_indices_vec.size() << "), "
                   "parameters_vec (" << parameters_vec.size() << ") is not matched "
               << "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        // Hinckley et al., J. Chem. Phys. (2013)
        std::string potential_formula ="energy;"
            "energy   = fdt3*fdtCS*attr/2;"  // Eq. (10)
            "attr     = epsilon*(1 - exp(-alpha_CS*dr))^2 * step(dr) - epsilon;"
            "fdt3     = max(f1*rect0t3,  rect1t3);"
            "fdtCS    = max(f2*rect0tCS, rect1tCS);"
            "rect0t3  = step(pi   + dt3)  * step(pi   - dt3);"
            "rect0tCS = step(pi   + dtCS) * step(pi   - dtCS);"
            "rect1t3  = step(pi/2 + dt3)  * step(pi/2 - dt3);"
            "rect1tCS = step(pi/2 + dtCS) * step(pi/2 - dtCS);"
            "f1       = 1 - cos(dt3)^2;"
            "f2       = 1 - cos(dtCS)^2;"
            "dr       = distance(d1, a3) - r0;"
            "dt3      = K_BP*(t3  - t03);"
            "dtCS     = K_CS*(tCS - t0CS);"
            "tCS      = angle(d2, d1, a3);"
            "t3       = acos(cost3lim);"
            "cost3lim = min(max(cost3, -0.99), 0.99);"
            "cost3    = sin(t1)*sin(t2)*cos(phi) - cos(t1)*cos(t2);"
            "t1       = angle(d2, d1, a1);"
            "t2       = angle(d1, a1, a2);"
            "phi      = dihedral(d2, d1, a1, a2);";

        const std::map<std::string, std::string> ff_params =
        {
            {"epsilon",  "TSPN2_CS" + ffgen_id_str_ + "_epsilon"},
            {"r0",       "TSPN2_CS" + ffgen_id_str_ + "_r0"},
            {"t03",      "TSPN2_CS" + ffgen_id_str_ + "_t03"},
            {"t0CS",     "TSPN2_CS" + ffgen_id_str_ + "_t0CS"},
            {"K_BP",     "TSPN2_CS" + ffgen_id_str_ + "_K_BP"},
            {"K_CS",     "TSPN2_CS" + ffgen_id_str_ + "_K_CS"},
            {"alpha_CS", "TSPN2_CS" + ffgen_id_str_ + "_alpha_CS"}
        };

        for (auto itr=ff_params.begin(); itr!=ff_params.end(); ++itr)
        {
            potential_formula = std::regex_replace(
                potential_formula, std::regex(itr->first), itr->second);
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

        if (0 < parameters_vec_.size())
        {
            chbond_ff->setCutoffDistance(parameters_vec_.at(0).cutoff());
        }

        chbond_ff->addPerAcceptorParameter(ff_params.at("epsilon"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("r0"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("t03"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("t0CS"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("K_BP"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("K_CS"));
        chbond_ff->addPerAcceptorParameter(ff_params.at("alpha_CS"));
        chbond_ff->addGlobalParameter("pi", Constant::pi);

        for (size_t i=0; i < donor_indices_vec_.size(); ++i)
        {
            const size_t d1_base_0 = donor_indices_vec_.at(i).at(0);
            const size_t d2_sugar  = donor_indices_vec_.at(i).at(1);
            const size_t d3_base_n = donor_indices_vec_.at(i).at(2);
            chbond_ff->addDonor(d1_base_0, d2_sugar, d3_base_n);
        }

        for (size_t i=0; i < acceptor_indices_vec_.size(); ++i)
        {
            const size_t a1_base_p = acceptor_indices_vec_.at(i).at(0);
            const size_t a2_sugar  = acceptor_indices_vec_.at(i).at(1);
            const size_t a3_base_c = acceptor_indices_vec_.at(i).at(2);

            const std::vector<double> parameters = {
                parameters_vec_.at(i).epsilon_CS(),
                parameters_vec_.at(i).r0_CS(),
                parameters_vec_.at(i).theta3_0(),
                parameters_vec_.at(i).theta_CS_0(),
                parameters_vec_.at(i).K_BP(),
                parameters_vec_.at(i).K_CS(),
                parameters_vec_.at(i).alpha_CS()
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

    std::string name() const noexcept {return "3SPN2CrossStacking";};

  private:
    std::vector<indices_type>   donor_indices_vec_;
    std::vector<indices_type>   acceptor_indices_vec_;
    index_pairs_type            ignore_list_;
    std::vector<parameter_list> parameters_vec_;
    const bool                  use_periodic_;
    const std::string           ffgen_id_str_;
};

// ----------------------------------------------------------------------------
// 3SPN2 parameter set
//
// The parameters are derived from the following literature
// - W. Lu, C. Bueno, N. P. Schafer, J. Moller, S. Jin, X. Chen, M. Chen, X. Gu,
//   A. Davtyan, J. J. de Pablo, and P. G. Wolynes, PLOS Comput. Biol. (2021)
//
//   - Table 11
//
// Also see
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// - G. S. Freeman, D. M. Hinckley, J. P. Lequieu, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2014)

// Parameter for 3SPN2
template<typename realT>
struct ThreeSPN2CrossStackingPotentialParameter
{
    using real_type = realT;
    real_type cutoff       = 18.0 * OpenMM::NmPerAngstrom; // [nm]
    real_type alpha_CS     =  4.0 / OpenMM::NmPerAngstrom;
    real_type K_CS         =  8.0;
    real_type K_BP         = 12.0;

    const std::map<std::string, real_type> theta3_0 = { // [radian]
        {"AT", 116.09 * OpenMM::RadiansPerDegree},
        {"TA", 116.09 * OpenMM::RadiansPerDegree},
        {"GC", 124.94 * OpenMM::RadiansPerDegree},
        {"CG", 124.94 * OpenMM::RadiansPerDegree}
    };

    const std::map<std::string, real_type> epsilon_CS = { //[kJ/mol]
        // sense strand
        {"AA5", 2.186}, {"AT5", 2.774}, {"AG5", 2.833}, {"AC5", 1.951},
        {"TA5", 2.774}, {"TT5", 2.186}, {"TG5", 2.539}, {"TC5", 2.980},
        {"GA5", 2.833}, {"GT5", 2.539}, {"GG5", 3.774}, {"GC5", 1.129},
        {"CA5", 1.951}, {"CT5", 2.980}, {"CG5", 1.129}, {"CC5", 4.802},
        // anti-sense strand
        {"AA3", 2.186}, {"AT3", 2.774}, {"AG3", 2.980}, {"AC3", 2.539},
        {"TA3", 2.774}, {"TT3", 2.186}, {"TG3", 1.951}, {"TC3", 2.833},
        {"GA3", 2.980}, {"GT3", 1.951}, {"GG3", 4.802}, {"GC3", 1.129},
        {"CA3", 2.539}, {"CT3", 2.833}, {"CG3", 1.129}, {"CC3", 3.774}
    };

    const std::map<std::string, real_type> r0_CS = { // [nm]
        // sense strand
        {"AA5", 6.208 * OpenMM::NmPerAngstrom},
        {"AT5", 6.876 * OpenMM::NmPerAngstrom},
        {"AG5", 6.072 * OpenMM::NmPerAngstrom},
        {"AC5", 6.811 * OpenMM::NmPerAngstrom},
        {"TA5", 6.876 * OpenMM::NmPerAngstrom},
        {"TT5", 7.480 * OpenMM::NmPerAngstrom},
        {"TG5", 6.771 * OpenMM::NmPerAngstrom},
        {"TC5", 7.453 * OpenMM::NmPerAngstrom},
        {"GA5", 6.072 * OpenMM::NmPerAngstrom},
        {"GT5", 6.771 * OpenMM::NmPerAngstrom},
        {"GG5", 5.921 * OpenMM::NmPerAngstrom},
        {"GC5", 6.688 * OpenMM::NmPerAngstrom},
        {"CA5", 6.811 * OpenMM::NmPerAngstrom},
        {"CT5", 7.453 * OpenMM::NmPerAngstrom},
        {"CG5", 6.688 * OpenMM::NmPerAngstrom},
        {"CC5", 7.409 * OpenMM::NmPerAngstrom},
        // anti-sense strand
        {"AA3", 5.435 * OpenMM::NmPerAngstrom},
        {"AT3", 6.295 * OpenMM::NmPerAngstrom},
        {"AG3", 5.183 * OpenMM::NmPerAngstrom},
        {"AC3", 6.082 * OpenMM::NmPerAngstrom},
        {"TA3", 6.295 * OpenMM::NmPerAngstrom},
        {"TT3", 7.195 * OpenMM::NmPerAngstrom},
        {"TG3", 6.028 * OpenMM::NmPerAngstrom},
        {"TC3", 6.981 * OpenMM::NmPerAngstrom},
        {"GA3", 5.183 * OpenMM::NmPerAngstrom},
        {"GT3", 6.028 * OpenMM::NmPerAngstrom},
        {"GG3", 4.934 * OpenMM::NmPerAngstrom},
        {"GC3", 5.811 * OpenMM::NmPerAngstrom},
        {"CA3", 6.082 * OpenMM::NmPerAngstrom},
        {"CT3", 6.981 * OpenMM::NmPerAngstrom},
        {"CG3", 5.811 * OpenMM::NmPerAngstrom},
        {"CC3", 6.757 * OpenMM::NmPerAngstrom}
    };

    const std::map<std::string, real_type> theta_CS_0 = { // [degree]
        // sense strand ("XY5" X=B0, Y=Bc)
        {"AA5", 154.38 * OpenMM::RadiansPerDegree},
        {"AT5", 159.10 * OpenMM::RadiansPerDegree},
        {"AG5", 152.46 * OpenMM::RadiansPerDegree},
        {"AC5", 158.38 * OpenMM::RadiansPerDegree},
        {"TA5", 147.10 * OpenMM::RadiansPerDegree},
        {"TT5", 153.79 * OpenMM::RadiansPerDegree},
        {"TG5", 144.44 * OpenMM::RadiansPerDegree},
        {"TC5", 151.48 * OpenMM::RadiansPerDegree},
        {"GA5", 154.69 * OpenMM::RadiansPerDegree},
        {"GT5", 157.83 * OpenMM::RadiansPerDegree},
        {"GG5", 153.43 * OpenMM::RadiansPerDegree},
        {"GC5", 158.04 * OpenMM::RadiansPerDegree},
        {"CA5", 152.99 * OpenMM::RadiansPerDegree},
        {"CT5", 159.08 * OpenMM::RadiansPerDegree},
        {"CG5", 150.53 * OpenMM::RadiansPerDegree},
        {"CC5", 157.17 * OpenMM::RadiansPerDegree},
        // anti-sense strand
        {"AA3", 116.88 * OpenMM::RadiansPerDegree},
        {"AT3", 121.74 * OpenMM::RadiansPerDegree},
        {"AG3", 114.23 * OpenMM::RadiansPerDegree},
        {"AC3", 119.06 * OpenMM::RadiansPerDegree},
        {"TA3", 109.42 * OpenMM::RadiansPerDegree},
        {"TT3", 112.95 * OpenMM::RadiansPerDegree},
        {"TG3", 107.32 * OpenMM::RadiansPerDegree},
        {"TC3", 110.56 * OpenMM::RadiansPerDegree},
        {"GA3", 119.34 * OpenMM::RadiansPerDegree},
        {"GT3", 124.72 * OpenMM::RadiansPerDegree},
        {"GG3", 116.51 * OpenMM::RadiansPerDegree},
        {"GC3", 121.98 * OpenMM::RadiansPerDegree},
        {"CA3", 114.60 * OpenMM::RadiansPerDegree},
        {"CT3", 118.26 * OpenMM::RadiansPerDegree},
        {"CG3", 112.45 * OpenMM::RadiansPerDegree},
        {"CC3", 115.88 * OpenMM::RadiansPerDegree}
    };
};

// parameter for 3SPN2C
template<typename realT>
struct ThreeSPN2CCrossStackingPotentialParameter
{
    using real_type = realT;
    real_type cutoff       = 18.0 * OpenMM::NmPerAngstrom; // [nm]
    real_type alpha_CS     =  4.0 / OpenMM::NmPerAngstrom;
    real_type K_CS         =  8.0;
    real_type K_BP         = 12.0;

    const std::map<std::string, real_type> theta3_0 = { // [radian]
        {"AT", 110.92 * OpenMM::RadiansPerDegree},
        {"TA", 110.92 * OpenMM::RadiansPerDegree},
        {"GC", 120.45 * OpenMM::RadiansPerDegree},
        {"CG", 120.45 * OpenMM::RadiansPerDegree},
    };

    const std::map<std::string, real_type> epsilon_CS = { // [kJ/mol]
        // sense strand
        {"AA5", 1.882},  {"AT5", 2.388},  {"AG5", 2.439},  {"AC5", 1.680},
        {"TA5", 2.388},  {"TT5", 1.882},  {"TG5", 2.187},  {"TC5", 2.566},
        {"GA5", 2.439},  {"GT5", 2.187},  {"GG5", 3.250},  {"GC5", 0.972},
        {"CA5", 1.680},  {"CT5", 2.566},  {"CG5", 0.972},  {"CC5", 4.135},
        // anti-sense strand
        {"AA3", 1.882},  {"AT3", 2.388},  {"AG3", 2.566},  {"AC3", 2.187},
        {"TA3", 2.388},  {"TT3", 1.882},  {"TG3", 1.680},  {"TC3", 2.439},
        {"GA3", 2.566},  {"GT3", 1.680},  {"GG3", 4.135},  {"GC3", 0.972},
        {"CA3", 2.187},  {"CT3", 2.439},  {"CG3", 0.972},  {"CC3", 3.250}
    };

    const std::map<std::string, real_type> r0_CS = { // [nm]
        // sense strand
        {"AA5", 6.420 * OpenMM::NmPerAngstrom},
        {"AT5", 6.770 * OpenMM::NmPerAngstrom},
        {"AG5", 6.270 * OpenMM::NmPerAngstrom},
        {"AC5", 6.840 * OpenMM::NmPerAngstrom},
        {"TA5", 6.770 * OpenMM::NmPerAngstrom},
        {"TT5", 7.210 * OpenMM::NmPerAngstrom},
        {"TG5", 6.530 * OpenMM::NmPerAngstrom},
        {"TC5", 7.080 * OpenMM::NmPerAngstrom},
        {"GA5", 6.270 * OpenMM::NmPerAngstrom},
        {"GT5", 6.530 * OpenMM::NmPerAngstrom},
        {"GG5", 5.740 * OpenMM::NmPerAngstrom},
        {"GC5", 6.860 * OpenMM::NmPerAngstrom},
        {"CA5", 6.840 * OpenMM::NmPerAngstrom},
        {"CT5", 7.080 * OpenMM::NmPerAngstrom},
        {"CG5", 6.860 * OpenMM::NmPerAngstrom},
        {"CC5", 6.790 * OpenMM::NmPerAngstrom},
        // anti-sense strand
        {"AA3", 5.580 * OpenMM::NmPerAngstrom},
        {"AT3", 6.140 * OpenMM::NmPerAngstrom},
        {"AG3", 5.630 * OpenMM::NmPerAngstrom},
        {"AC3", 6.180 * OpenMM::NmPerAngstrom},
        {"TA3", 6.140 * OpenMM::NmPerAngstrom},
        {"TT3", 6.800 * OpenMM::NmPerAngstrom},
        {"TG3", 6.070 * OpenMM::NmPerAngstrom},
        {"TC3", 6.640 * OpenMM::NmPerAngstrom},
        {"GA3", 5.630 * OpenMM::NmPerAngstrom},
        {"GT3", 6.070 * OpenMM::NmPerAngstrom},
        {"GG3", 5.870 * OpenMM::NmPerAngstrom},
        {"GC3", 5.660 * OpenMM::NmPerAngstrom},
        {"CA3", 6.180 * OpenMM::NmPerAngstrom},
        {"CT3", 6.640 * OpenMM::NmPerAngstrom},
        {"CG3", 5.660 * OpenMM::NmPerAngstrom},
        {"CC3", 6.800 * OpenMM::NmPerAngstrom}
    };

    const std::map<std::string, real_type> theta_CS_0 = {
        // sense strand ("XY5", X=B0, Y=Bc)
        {"AA5", 154.04 * OpenMM::RadiansPerDegree},
        {"AT5", 158.77 * OpenMM::RadiansPerDegree},
        {"AG5", 153.88 * OpenMM::RadiansPerDegree},
        {"AC5", 157.69 * OpenMM::RadiansPerDegree},
        {"TA5", 148.62 * OpenMM::RadiansPerDegree},
        {"TT5", 155.05 * OpenMM::RadiansPerDegree},
        {"TG5", 147.54 * OpenMM::RadiansPerDegree},
        {"TC5", 153.61 * OpenMM::RadiansPerDegree},
        {"GA5", 153.91 * OpenMM::RadiansPerDegree},
        {"GT5", 155.72 * OpenMM::RadiansPerDegree},
        {"GG5", 151.84 * OpenMM::RadiansPerDegree},
        {"GC5", 157.80 * OpenMM::RadiansPerDegree},
        {"CA5", 152.04 * OpenMM::RadiansPerDegree},
        {"CT5", 157.72 * OpenMM::RadiansPerDegree},
        {"CG5", 151.65 * OpenMM::RadiansPerDegree},
        {"CC5", 154.49 * OpenMM::RadiansPerDegree},
        // anti-sense strand
        {"AA3", 116.34 * OpenMM::RadiansPerDegree},
        {"AT3", 119.61 * OpenMM::RadiansPerDegree},
        {"AG3", 115.19 * OpenMM::RadiansPerDegree},
        {"AC3", 120.92 * OpenMM::RadiansPerDegree},
        {"TA3", 107.40 * OpenMM::RadiansPerDegree},
        {"TT3", 110.76 * OpenMM::RadiansPerDegree},
        {"TG3", 106.33 * OpenMM::RadiansPerDegree},
        {"TC3", 111.57 * OpenMM::RadiansPerDegree},
        {"GA3", 121.61 * OpenMM::RadiansPerDegree},
        {"GT3", 124.92 * OpenMM::RadiansPerDegree},
        {"GG3", 120.52 * OpenMM::RadiansPerDegree},
        {"GC3", 124.88 * OpenMM::RadiansPerDegree},
        {"CA3", 112.45 * OpenMM::RadiansPerDegree},
        {"CT3", 115.43 * OpenMM::RadiansPerDegree},
        {"CG3", 110.51 * OpenMM::RadiansPerDegree},
        {"CC3", 115.80 * OpenMM::RadiansPerDegree}
    };
};

#endif // OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP
