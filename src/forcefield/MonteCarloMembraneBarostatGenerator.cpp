#include "MonteCarloMembraneBarostatGenerator.hpp"

#include <iostream>
#include <iomanip>

std::unique_ptr<OpenMM::Force> MonteCarloMembraneBarostatGenerator::generate() const
{
    auto barostat =
        std::make_unique<OpenMM::MonteCarloMembraneBarostat>(
            default_pressure_, default_surface_tension_, temperature_,
            xymode_, zmode_);
    return barostat;
}
void MonteCarloMembraneBarostatGenerator::dump_info() const
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

