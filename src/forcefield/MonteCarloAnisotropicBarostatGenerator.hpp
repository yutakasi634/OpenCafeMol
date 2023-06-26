#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP

#include <OpenMM.h>

class MonteCarloAnisotropicBarostatGenerator
{
  public:
    MonteCarloAnisotropicBarostatGenerator(
        const std::array<bool, 3>   scale_axis,       const double temperature,
        const std::array<double, 3> default_pressure, const std::size_t frequency)
        : scale_axis_(scale_axis), temperature_(temperature), default_pressure_(default_pressure),
          frequency_(frequency)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept
    {
        auto barostat =
            std::make_unique<OpenMM::MonteCarloAnisotropicBarostat>(
                OpenMM::Vec3(default_pressure_[0], default_pressure_[1], default_pressure_[2]),
                temperature_, scale_axis_[0], scale_axis_[1], scale_axis_[2], frequency_);

        return barostat;
    }

    const std::array<bool, 3>&   scale_axis()       const noexcept { return scale_axis_; }
    const std::array<double, 3>& default_pressure() const noexcept { return default_pressure_; }

    double      temperature() const noexcept { return temperature_; }
    std::size_t frequency()   const noexcept { return frequency_; }

  private:
    std::array<bool, 3>   scale_axis_;
    double                temperature_;
    std::array<double, 3> default_pressure_;
    std::size_t           frequency_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
