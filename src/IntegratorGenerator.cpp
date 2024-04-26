#include "IntegratorGenerator.hpp"

#include <iomanip>
#include <sstream>

std::unique_ptr<OpenMM::Integrator> LangevinIntegratorGenerator::generate() const
{
    auto integrator = std::make_unique<OpenMM::LangevinIntegrator>(
            this->temperature_,
            this->friction_coeff_,
            this->delta_t_);
    integrator->setRandomNumberSeed(this->seed_);

    return integrator;
}

std::string LangevinIntegratorGenerator::dump_info() const
{
    std::ostringstream oss;

    oss << "    temperature : "
        << std::setw(7) << std::fixed << std::setprecision(2)
        << this->temperature_ << " K\n";
    oss << "    delta t     : "
        << std::setw(7) << std::fixed << std::setprecision(3)
        << this->delta_t_ << " ps\n";
    oss << "    gamma       : "
        << std::setw(7) << std::fixed << std::setprecision(3)
        << this->friction_coeff_ << " ps^-1\n";

    if(this->seed_ == 0)
    {
        oss << "    seed        : "
            << "not specified or 0. random seed will be chosen.\n";
    }
    else
    {
        oss << "    seed        : "
            << std::setw(7) << this->seed_ << '\n';
    }
    return oss.str();
}


