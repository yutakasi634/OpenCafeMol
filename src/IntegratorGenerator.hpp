#ifndef OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP
#define OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP

#include <OpenMM.h>

#include <memory>
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

    std::unique_ptr<OpenMM::Integrator> generate() const override;

    std::string dump_info() const override;

    double temperature() const override { return temperature_; }

  private:

    double temperature_;
    double friction_coeff_;
    double delta_t_;
    int    seed_;
};

#endif// OPEN_AICG2_PLUS_INTEGRATOR_GENERATOR_HPP
