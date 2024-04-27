#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP

#include "BarostatGeneratorBase.hpp"

#include <OpenMM.h>

#include <array>

class MonteCarloAnisotropicBarostatGenerator final : public BarostatGeneratorBase
{
  public:
    MonteCarloAnisotropicBarostatGenerator(
        const std::array<bool, 3>   scale_axis,       const double temperature,
        const std::array<double, 3> default_pressure, const std::size_t frequency)
        : scale_axis_(scale_axis), temperature_(temperature), default_pressure_(default_pressure),
          frequency_(frequency)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    const std::array<bool, 3>&   scale_axis()       const noexcept { return scale_axis_; }
    const std::array<double, 3>& default_pressure() const noexcept { return default_pressure_; }

    double      temperature() const noexcept override { return temperature_; }
    std::size_t frequency()   const noexcept override { return frequency_; }
    std::string name()        const noexcept override { return "MonteCarloAnisotropicBarostat"; }

    void dump_info() const noexcept override;

  private:
    std::array<bool, 3>   scale_axis_;
    double                temperature_;
    std::array<double, 3> default_pressure_;
    std::size_t           frequency_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
