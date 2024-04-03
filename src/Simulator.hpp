#ifndef OPEN_AICG2_PLUS_SIMULATOR_HPP
#define OPEN_AICG2_PLUS_SIMULATOR_HPP

#include <OpenMM.h>
#include "util/Constants.hpp"
#include "Observer.hpp"
#include "IntegratorGenerator.hpp"

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
        bool energy_minimization = false, bool output_progress = true)
        : system_ptr_(system_gen.generate().release()),
          integrator_gen_ptr_(std::move(integrator_gen_ptr)),
          integrator_ptr_(integrator_gen_ptr_->generate().release()),
          context_(*system_ptr_, *integrator_ptr_, platform, platform_properties),
          total_step_(total_step), save_step_(save_step),
          energy_minimization_(energy_minimization), output_progress_(output_progress),
          progress_bar_(/* width of bar = */ 50)
    {
        std::cerr << "initializing simulator..." << std::endl;
        std::cerr << "    total step : " << total_step_ << " step" << std::endl;
        std::cerr << "    save step  : " << save_step_  << " step" << std::endl;
        std::cerr << "integrator information" << std::endl;
        std::cerr << integrator_gen_ptr_->dump_info();

        std::cerr << "initializing observers..." << std::endl;
        for(auto& observer : observers)
        {
            observers_.push_back(std::move(observer));
            observers_.back()->initialize(system_ptr_);
        }
    }

    void initialize(const std::vector<OpenMM::Vec3>& initial_position)
    {
        // Set starting positions of the atoms.
        context_.setPositions(initial_position);
        if(energy_minimization_)
        {
            std::cerr << "energy minimization in progress..." << std::endl;
            OpenMM::LocalEnergyMinimizer::minimize(context_);
        }
    }

    void set_velocity()
    {
        context_.setVelocitiesToTemperature(integrator_gen_ptr_->temperature());
    }

    void set_velocity(const std::vector<OpenMM::Vec3>& initial_velocity)
    {
        context_.setVelocities(initial_velocity);
    }

    void run()
    {
        const std::size_t total_frame = std::floor(total_step_/save_step_);
        for(std::size_t frame_num=0; frame_num<total_frame; ++frame_num)
        {
            const std::size_t step = frame_num*save_step_;
            for(const auto& observer : observers_)
            {
                observer->output(frame_num*save_step_, context_);
            }

            integrator_ptr_->step(save_step_);

            if(output_progress_)
            {
                progress_bar_.format(step, total_step_, std::cerr);
            }
        }

        for(const auto& observer : observers_)
        {
            observer->finalize();
        }
        if(output_progress_)
        {
            progress_bar_.finalize(std::cerr);
        }
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
