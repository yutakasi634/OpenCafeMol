#ifndef OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP
#define OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP

#include <OpenMM.h>

#include <string>
#include <map>

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

#endif// OPEN_AICG2_PLUS_3SPN2_DEFAULT_PARAMETERS_HPP
