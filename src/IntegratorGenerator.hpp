#ifndef OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP
#define OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP

#include <OpenMM.h>

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

class IntegratorGeneratorBase
{
  public:
    virtual ~IntegratorGeneratorBase() = default;

    virtual std::unique_ptr<OpenMM::Integrator> generate() const = 0;
    virtual std::string dump_info() const = 0;
    virtual double temperature() const = 0;
};

class LangevinIntegratorGenerator final: public IntegratorGeneratorBase
{
  public:

    LangevinIntegratorGenerator(const double temperature,
        const double friction_coeff, const double delta_t, const int seed)
        : temperature_(temperature), friction_coeff_(friction_coeff),
          delta_t_(delta_t), seed_(seed)
    {}
    ~LangevinIntegratorGenerator() override = default;

    std::unique_ptr<OpenMM::Integrator> generate() const override
    {
        auto integrator = std::make_unique<OpenMM::LangevinIntegrator>(
                this->temperature_,
                this->friction_coeff_,
                this->delta_t_);
        integrator->setRandomNumberSeed(this->seed_);

        return integrator;
    }

    std::string dump_info() const override
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

    double temperature() const override { return temperature_; }

  private:

    double temperature_;
    double friction_coeff_;
    double delta_t_;
    int    seed_;
};

#endif// OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP
