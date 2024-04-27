#include "IntegratorGenerator.hpp"

#include <fmt/core.h>

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
    std::string info;
    info += fmt::format("    temperature : {:7.2f} K\n",     this->temperature_);
    info += fmt::format("    delta_t     : {:7.3f} ps\n",    this->delta_t_);
    info += fmt::format("    gamma       : {:7.3f} ps^-1\n", this->friction_coeff_);


    std::string seed;
    if(this->seed_ == 0)
    {
        seed = "not specified or 0. random seed will be chosen.";
    }
    else
    {
        seed = fmt::format("{:7}", this->seed_);
    }
    info += fmt::format("    seed        : {}\n", seed);
    return info;
}


