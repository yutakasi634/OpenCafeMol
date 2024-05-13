#ifndef OPEN_AICG2_PLUS_SIMULATOR_HPP
#define OPEN_AICG2_PLUS_SIMULATOR_HPP

#include "util/ProgressBar.hpp"
#include "Observer.hpp"
#include "IntegratorGenerator.hpp"

#include <OpenMM.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

class Simulator
{
  public:
    Simulator(
        const SystemGenerator& system_gen,
        std::unique_ptr<IntegratorGeneratorBase>&& integrator_gen_ptr,
        OpenMM::Platform& platform,
        const std::map<std::string, std::string>& platform_properties,
        const std::vector<OpenMM::Vec3>& initial_position,
        const std::size_t total_step, const std::size_t save_step,
        std::vector<std::unique_ptr<ObserverBase>>& observers,
        bool energy_minimization = false, bool output_progress = true)
        : Simulator(system_gen, std::move(integrator_gen_ptr),
                    platform, platform_properties,
                    total_step, save_step, observers,
                    energy_minimization, output_progress)
    {
        this->initialize(initial_position);
    }

    Simulator(
        const SystemGenerator& system_gen,
        std::unique_ptr<IntegratorGeneratorBase>&& integrator_gen_ptr,
        OpenMM::Platform& platform,
        const std::map<std::string, std::string>& platform_properties,
        const std::size_t total_step, const std::size_t save_step,
        std::vector<std::unique_ptr<ObserverBase>>& observers,
        bool energy_minimization = false,
        bool output_progress = true);

    void initialize(const std::vector<OpenMM::Vec3>& initial_position);
    void run();

    void set_velocity(const std::int64_t vel_seed)
    {
        if(vel_seed == 0)
        {
            context_.setVelocitiesToTemperature(integrator_gen_ptr_->temperature());
        }
        else
        {
            context_.setVelocitiesToTemperature(integrator_gen_ptr_->temperature(),
                                                static_cast<int>(vel_seed));
        }
    }

    void set_velocity(const std::vector<OpenMM::Vec3>& initial_velocity)
    {
        context_.setVelocities(initial_velocity);
    }

    std::size_t total_step() const noexcept { return total_step_; }
    std::size_t save_step()  const noexcept { return save_step_; }

  private:
    std::unique_ptr<OpenMM::System>            system_ptr_;
    std::unique_ptr<IntegratorGeneratorBase>   integrator_gen_ptr_;
    std::unique_ptr<OpenMM::Integrator>        integrator_ptr_;
    OpenMM::Context                            context_;
    std::size_t                                total_step_;
    std::size_t                                save_step_;
    std::vector<std::unique_ptr<ObserverBase>> observers_;
    bool                                       energy_minimization_;
    bool                                       output_progress_;
    ProgressBar                                progress_bar_;
};

#endif // OPEN_AICG2_PLUS_SIMULATOR_HPP
