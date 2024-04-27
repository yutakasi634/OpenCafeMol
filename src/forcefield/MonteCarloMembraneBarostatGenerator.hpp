#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP

#include "BarostatGeneratorBase.hpp"

#include <OpenMM.h>

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

    std::unique_ptr<OpenMM::Force> generate() const override;

    double      temperature() const override { return temperature_; }
    std::size_t frequency()   const override { return frequency_; }
    std::string name()        const override { return "MonteCarloMembraneBarostat"; }

    void dump_info() const override;

  private:
    double      default_pressure_;        // bar
    double      default_surface_tension_; // bar*nm
    double      temperature_;             // K
    xymode_type xymode_;
    zmode_type  zmode_;
    std::size_t frequency_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_MEMBRANE_BAROSTAT_GENERATOR_HPP
