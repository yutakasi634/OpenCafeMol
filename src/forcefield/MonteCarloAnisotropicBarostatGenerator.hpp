#ifndef OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
#define OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP

#include "BarostatGeneratorBase.hpp"

#include <OpenMM.h>

#include <array>
#include <iostream>
#include <iomanip>
#include <sstream>

class MonteCarloAnisotropicBarostatGenerator final : public BarostatGeneratorBase
{
  public:
    MonteCarloAnisotropicBarostatGenerator(
        const std::array<bool, 3>   scale_axis,       const double temperature,
        const std::array<double, 3> default_pressure, const std::size_t frequency)
        : scale_axis_(scale_axis), temperature_(temperature), default_pressure_(default_pressure),
          frequency_(frequency)
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        auto barostat =
            std::make_unique<OpenMM::MonteCarloAnisotropicBarostat>(
                OpenMM::Vec3(default_pressure_[0], default_pressure_[1], default_pressure_[2]),
                temperature_, scale_axis_[0], scale_axis_[1], scale_axis_[2], frequency_);

        return barostat;
    }

    const std::array<bool, 3>&   scale_axis()       const noexcept { return scale_axis_; }
    const std::array<double, 3>& default_pressure() const noexcept { return default_pressure_; }

    double      temperature() const noexcept override { return temperature_; }
    std::size_t frequency()   const noexcept override { return frequency_; }
    std::string name()        const noexcept override { return "MonteCarloAnisotropicBarostat"; }

    void dump_info() const noexcept override
    {
        std::cerr << "        scaling axis              : ";
        std::stringstream scaling_axis;
        if(scale_axis_[0]){ scaling_axis << "X"; }
        if(scale_axis_[1]){ scaling_axis << "Y"; }
        if(scale_axis_[2]){ scaling_axis << "Z"; }
        std::cerr << std::setw(14) << scaling_axis.str() << std::endl;

        std::cerr << "        default pressure"
            << std::fixed << std::setprecision(2) << std::endl;
        if(scale_axis_[0]){ std::cerr << "            X: "
            << std::setw(7) << default_pressure_[0] << std::endl; }
        if(scale_axis_[1]){ std::cerr << "            Y: "
            << std::setw(7) << default_pressure_[1] << std::endl; }
        if(scale_axis_[2]){ std::cerr << "            Z: "
            << std::setw(7) << default_pressure_[2] << std::endl; }
    }

  private:
    std::array<bool, 3>   scale_axis_;
    double                temperature_;
    std::array<double, 3> default_pressure_;
    std::size_t           frequency_;
};

#endif // OPEN_AICG2_PLUS_MONTE_CARLO_ANISOTROPIC_BAROSTAT_GENERATOR_HPP
