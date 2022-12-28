#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP

class MonteCarloAnisotropicBarostatGenerator
{
  public:
    MonteCarloAnisotropicBarostatGenerator(
        const std::array<bool, 3>   scale_axis, const double temperature,
        const std::array<double, 3> default_pressure)
        : scale_axis_(scale_axis), temperature_(temperature), default_pressure_(default_pressure)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept
    {
        auto barostat =
            std::make_unique<OpenMM::MonteCarloAnisotropicBarostat>(
                    OpenMM::Vec3(default_pressure_[0], default_pressure_[1], default_pressure_[2]),
                    temperature_, scale_axis_[0], scale_axis_[1], scale_axis_[2]);

        return barostat;
    }

  private:
    const std::array<bool, 3>   scale_axis_;
    const double                temperature_;
    const std::array<double, 3> default_pressure_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
