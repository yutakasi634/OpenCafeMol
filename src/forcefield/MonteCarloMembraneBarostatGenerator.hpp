#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP

#include "BarostatGeneratorBase.hpp"

#include <OpenMM.h>

#include <iostream>
#include <iomanip>
#include <memory>
#include <string>

class MonteCarloMembraneBarostatGenerator final : public BarostatGeneratorBase
{
  public:
    using xymode_type = OpenMM::MonteCarloMembraneBarostat::XYMode;
    using zmode_type  = OpenMM::MonteCarloMembraneBarostat::ZMode;

  public:
    MonteCarloMembraneBarostatGenerator(
        const double default_pressure, const double default_surface_tension,
        const double temperature,      const xymode_type xymode, const zmode_type zmode,
        const std::size_t frequency)
        : default_pressure_(default_pressure),
          default_surface_tension_(default_surface_tension), temperature_(temperature),
          xymode_(xymode), zmode_(zmode), frequency_(frequency)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        auto barostat =
            std::make_unique<OpenMM::MonteCarloMembraneBarostat>(
                default_pressure_, default_surface_tension_, temperature_,
                xymode_, zmode_);
        return barostat;
    }

    double      temperature() const noexcept override { return temperature_; }
    std::size_t frequency()   const noexcept override { return frequency_; }
    std::string name()        const noexcept override { return "MonteCarloMembraneBarostat"; }

    void dump_info() const noexcept override
    {
        std::cerr << "        XYMode                    : " << std::setw(14);
        if(xymode_ == xymode_type::XYIsotropic)
        {
            std::cerr << "Isotropic";
        }
        else if(xymode_ == xymode_type::XYAnisotropic)
        {
            std::cerr << "Anisotropic";
        }
        std::cerr << std::endl;

        std::cerr << "        ZMode                     : " << std::setw(14);
        if(zmode_ == zmode_type::ZFree)
        {
            std::cerr << "Free";
        }
        else if(zmode_ == zmode_type::ZFixed)
        {
            std::cerr << "Fixed";
        }
        else if(zmode_ == zmode_type::ConstantVolume)
        {
            std::cerr << "ConstantVolume";
        }
        std::cerr << std::endl;

        std::cerr << std::fixed << std::setprecision(2);
        std::cerr << "        default pressure          : "
            << std::setw(14) << default_pressure_ << " bar" << std::endl;
        std::cerr << "        default surface tention   : "
            << std::setw(14) << default_surface_tension_ << " bar nm" << std::endl;
    }

  private:
    double      default_pressure_;        // bar
    double      default_surface_tension_; // bar*nm
    double      temperature_;             // K
    xymode_type xymode_;
    zmode_type  zmode_;
    std::size_t frequency_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP
