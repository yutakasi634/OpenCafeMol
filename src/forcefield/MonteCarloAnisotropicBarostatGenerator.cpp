#include "MonteCarloAnisotropicBarostatGenerator.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>

std::unique_ptr<OpenMM::Force> MonteCarloAnisotropicBarostatGenerator::generate() const
{
    auto barostat =
        std::make_unique<OpenMM::MonteCarloAnisotropicBarostat>(
            OpenMM::Vec3(default_pressure_[0], default_pressure_[1], default_pressure_[2]),
            temperature_, scale_axis_[0], scale_axis_[1], scale_axis_[2], frequency_);

    return barostat;
}

void MonteCarloAnisotropicBarostatGenerator::dump_info() const
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


