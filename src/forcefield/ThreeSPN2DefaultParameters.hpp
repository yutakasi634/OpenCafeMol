#ifndef OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP
#define OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP

#include <OpenMM.h>

#include <string>
#include <map>

// =============================================================================
// base pair

// default parameter table
struct ThreeSPN2BasePairPotentialDefaultParameter
{
    virtual ~ThreeSPN2BasePairPotentialDefaultParameter() = default;

    virtual std::string name() const = 0;

    virtual double cutoff  () const = 0;
    virtual double alpha_BP() const = 0;
    virtual double K_BP    () const = 0;

    virtual std::map<std::string, double> epsilon_BP() const = 0;
    virtual std::map<std::string, double> r0()         const = 0;
    virtual std::map<std::string, double> theta0_1()   const = 0;
    virtual std::map<std::string, double> theta0_2()   const = 0;
    virtual std::map<std::string, double> phi0()       const = 0;
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

struct ThreeSPN2BasePairPotentialParameter final
    : public ThreeSPN2BasePairPotentialDefaultParameter
{
    ~ThreeSPN2BasePairPotentialParameter() override = default;

    std::string name    () const override { return "3SPN2";}
    double      cutoff  () const override { return 18.0 * OpenMM::NmPerAngstrom;} // [nm]
    double      alpha_BP() const override { return  2.0 / OpenMM::NmPerAngstrom;} // [1/nm]
    double      K_BP    () const override { return 12.0;}

    std::map<std::string, double> epsilon_BP() const override
    {
        return std::map<std::string, double>{ // [kJ/mol]
            {"AT", 16.73},
            {"TA", 16.73},
            {"GC", 21.18},
            {"CG", 21.18}
        };
    }
    std::map<std::string, double> r0() const override
    {
        return std::map<std::string, double>{ // [nm]
            {"AT", 5.941 * OpenMM::NmPerAngstrom},
            {"TA", 5.941 * OpenMM::NmPerAngstrom},
            {"GC", 5.530 * OpenMM::NmPerAngstrom},
            {"GC", 5.530 * OpenMM::NmPerAngstrom}
        };
    }
    std::map<std::string, double> theta0_1() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 156.54 * OpenMM::RadiansPerDegree},
            {"TA", 135.78 * OpenMM::RadiansPerDegree},
            {"GC", 159.81 * OpenMM::RadiansPerDegree},
            {"GC", 141.16 * OpenMM::RadiansPerDegree}
        };
    }
    std::map<std::string, double> theta0_2() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 135.78 * OpenMM::RadiansPerDegree},
            {"TA", 156.54 * OpenMM::RadiansPerDegree},
            {"GC", 141.16 * OpenMM::RadiansPerDegree},
            {"GC", 159.81 * OpenMM::RadiansPerDegree}
        };
    }
    std::map<std::string, double> phi0() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", -38.35 * OpenMM::RadiansPerDegree},
            {"TA", -38.35 * OpenMM::RadiansPerDegree},
            {"GC", -42.98 * OpenMM::RadiansPerDegree},
            {"GC", -42.98 * OpenMM::RadiansPerDegree}
        };
    }
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

struct ThreeSPN2CBasePairPotentialParameter final
    : public ThreeSPN2BasePairPotentialDefaultParameter
{
    ~ThreeSPN2CBasePairPotentialParameter() override = default;

    std::string name    () const override {return "3SPN2C";}
    double      cutoff  () const override {return 18.0 * OpenMM::NmPerAngstrom;}  // [nm]
    double      alpha_BP() const override {return  2.0 / OpenMM::NmPerAngstrom;}  // [1/nm]
    double      K_BP    () const override {return 12.0;}

    std::map<std::string, double> epsilon_BP() const override
    {
        return std::map<std::string, double>{ // [kJ/mol]
            {"AT", 14.41},
            {"TA", 14.41},
            {"GC", 18.24},
            {"CG", 18.24}
        };
    }
    std::map<std::string, double> r0() const override
    {
        return std::map<std::string, double>{ // [nm]
            {"AT", 5.82 * OpenMM::NmPerAngstrom},
            {"TA", 5.82 * OpenMM::NmPerAngstrom},
            {"GC", 5.52 * OpenMM::NmPerAngstrom},
            {"GC", 5.52 * OpenMM::NmPerAngstrom}
        };
    }
    std::map<std::string, double> theta0_1() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 153.17 * OpenMM::RadiansPerDegree},
            {"TA", 133.51 * OpenMM::RadiansPerDegree},
            {"GC", 159.50 * OpenMM::RadiansPerDegree},
            {"GC", 138.08 * OpenMM::RadiansPerDegree}
        };
    }
    std::map<std::string, double> theta0_2() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 133.51 * OpenMM::RadiansPerDegree},
            {"TA", 153.17 * OpenMM::RadiansPerDegree},
            {"GC", 138.08 * OpenMM::RadiansPerDegree},
            {"GC", 159.50 * OpenMM::RadiansPerDegree}
        };
    }
    std::map<std::string, double> phi0() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", -38.18 * OpenMM::RadiansPerDegree},
            {"TA", -38.18 * OpenMM::RadiansPerDegree},
            {"GC", -35.75 * OpenMM::RadiansPerDegree},
            {"GC", -35.75 * OpenMM::RadiansPerDegree}
        };
    }
};

// =============================================================================
// cross stacking


struct ThreeSPN2CrossStackingPotentialDefaultParameter
{
    virtual ~ThreeSPN2CrossStackingPotentialDefaultParameter() = default;

    virtual std::string name    () const = 0;
    virtual double      cutoff  () const = 0;
    virtual double      alpha_CS() const = 0;
    virtual double      K_CS    () const = 0;
    virtual double      K_BP    () const = 0;

    virtual std::map<std::string, double> theta3_0  () const = 0; // [radian]
    virtual std::map<std::string, double> epsilon_CS() const = 0; // [kJ/mol]
    virtual std::map<std::string, double> r0_CS     () const = 0; // [nm]
    virtual std::map<std::string, double> theta_CS_0() const = 0; // [degree]
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
struct ThreeSPN2CrossStackingPotentialParameter final
    : public ThreeSPN2CrossStackingPotentialDefaultParameter
{
    ~ThreeSPN2CrossStackingPotentialParameter() override = default;

    std::string name    () const override {return "3SPN2";}
    double      cutoff  () const override {return 18.0 * OpenMM::NmPerAngstrom;}
    double      alpha_CS() const override {return  4.0 / OpenMM::NmPerAngstrom;}
    double      K_CS    () const override {return  8.0;}
    double      K_BP    () const override {return 12.0;}

    std::map<std::string, double> theta3_0() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 116.09 * OpenMM::RadiansPerDegree},
            {"TA", 116.09 * OpenMM::RadiansPerDegree},
            {"GC", 124.94 * OpenMM::RadiansPerDegree},
            {"CG", 124.94 * OpenMM::RadiansPerDegree}
        };
    }

    std::map<std::string, double> epsilon_CS() const override
    {
        return std::map<std::string, double>{ //[kJ/mol]
            // sense strand ("XY5" X=B0, Y=Bc)
            {"AA5", 2.186}, {"AT5", 2.774}, {"AG5", 2.833}, {"AC5", 1.951},
            {"TA5", 2.774}, {"TT5", 2.186}, {"TG5", 2.539}, {"TC5", 2.980},
            {"GA5", 2.833}, {"GT5", 2.539}, {"GG5", 3.774}, {"GC5", 1.129},
            {"CA5", 1.951}, {"CT5", 2.980}, {"CG5", 1.129}, {"CC5", 4.802},
            // anti-sense strand ("XY5" X=B0, Y=Bc)
            {"AA3", 2.186}, {"AT3", 2.774}, {"AG3", 2.980}, {"AC3", 2.539},
            {"TA3", 2.774}, {"TT3", 2.186}, {"TG3", 1.951}, {"TC3", 2.833},
            {"GA3", 2.980}, {"GT3", 1.951}, {"GG3", 4.802}, {"GC3", 1.129},
            {"CA3", 2.539}, {"CT3", 2.833}, {"CG3", 1.129}, {"CC3", 3.774}
        };
    }

    std::map<std::string, double> r0_CS() const override
    {
        return std::map<std::string, double>{ // [nm]
            // sense strand ("XY5" X=B0, Y=Bc)
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
            // anti-sense strand ("XY5" X=B0, Y=Bc)
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
    }

    std::map<std::string, double> theta_CS_0() const override
    {
        return std::map<std::string, double>{ // [degree]
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
            // anti-sense strand ("XY5" X=B0, Y=Bc)
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
    }
};

// parameter for 3SPN2C
struct ThreeSPN2CCrossStackingPotentialParameter final
    : public ThreeSPN2CrossStackingPotentialDefaultParameter
{
    ~ThreeSPN2CCrossStackingPotentialParameter() override = default;

    std::string name    () const override {return "3SPN2C";}
    double      cutoff  () const override {return 18.0 * OpenMM::NmPerAngstrom;} // [nm]
    double      alpha_CS() const override {return  4.0 / OpenMM::NmPerAngstrom;}
    double      K_CS    () const override {return  8.0;}
    double      K_BP    () const override {return 12.0;}

    std::map<std::string, double> theta3_0() const override
    {
        return std::map<std::string, double>{ // [radian]
            {"AT", 110.92 * OpenMM::RadiansPerDegree},
            {"TA", 110.92 * OpenMM::RadiansPerDegree},
            {"GC", 120.45 * OpenMM::RadiansPerDegree},
            {"CG", 120.45 * OpenMM::RadiansPerDegree},
        };
    }
    std::map<std::string, double> epsilon_CS() const override
    {
        return std::map<std::string, double>{ // [kJ/mol]
            // sense strand ("XY5" X=B0, Y=Bc)
            {"AA5", 1.882},  {"AT5", 2.388},  {"AG5", 2.439},  {"AC5", 1.680},
            {"TA5", 2.388},  {"TT5", 1.882},  {"TG5", 2.187},  {"TC5", 2.566},
            {"GA5", 2.439},  {"GT5", 2.187},  {"GG5", 3.250},  {"GC5", 0.972},
            {"CA5", 1.680},  {"CT5", 2.566},  {"CG5", 0.972},  {"CC5", 4.135},
            // anti-sense strand ("XY5" X=B0, Y=Bc)
            {"AA3", 1.882},  {"AT3", 2.388},  {"AG3", 2.566},  {"AC3", 2.187},
            {"TA3", 2.388},  {"TT3", 1.882},  {"TG3", 1.680},  {"TC3", 2.439},
            {"GA3", 2.566},  {"GT3", 1.680},  {"GG3", 4.135},  {"GC3", 0.972},
            {"CA3", 2.187},  {"CT3", 2.439},  {"CG3", 0.972},  {"CC3", 3.250}
        };
    }
    std::map<std::string, double> r0_CS() const override
    {
        return std::map<std::string, double>{ // [nm]
            // sense strand ("XY5" X=B0, Y=Bc)
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
            // anti-sense strand ("XY5" X=B0, Y=Bc)
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
    }
    std::map<std::string, double> theta_CS_0() const override
    {
        return std::map<std::string, double>{
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
            // anti-sense strand ("XY5" X=B0, Y=Bc)
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
    }
};
#endif// OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP
