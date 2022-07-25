#ifndef OPEN_AICG2_PLUS_SIMULATOR_HPP
#define OPEN_AICG2_PLUS_SIMULATOR_HPP

#include <OpenMM.h>
#include "util/Constants.hpp"
#include "Observer.hpp"

class Simulator
{
  public:
    Simulator(
        std::unique_ptr<OpenMM::System>&& system_ptr, OpenMM::LangevinIntegrator&& integrator,
        const std::vector<OpenMM::Vec3>& initial_position,
        const std::size_t total_step, const std::size_t save_step,
        const Observer&& observer)
        : Simulator(std::move(system_ptr), std::move(integrator),
                    total_step, save_step, std::move(observer))
    {
        this->initialize(initial_position);
    }

    Simulator(
        std::unique_ptr<OpenMM::System>&& system_ptr, OpenMM::LangevinIntegrator&& integrator,
        const std::size_t total_step, const std::size_t save_step,
        const Observer&& observer)
        : system_ptr_(std::move(system_ptr)), integrator_(std::move(integrator)),
          context_(*system_ptr_, integrator_,
                    OpenMM::Platform::getPlatformByName("CUDA")),
          total_step_(total_step), save_step_(save_step),
          observer_(std::move(observer))
    {
        std::cerr << "initializing simulator..." << std::endl;
        std::cerr << "    total step : " << total_step_ << " step" << std::endl;
        std::cerr << "    save step  : " << save_step_  << " step" << std::endl;
        std::cerr << "integrator information" << std::endl;
        std::cerr << "    temperature : "
            << std::setw(7) << std::fixed << std::setprecision(2)
            << integrator_.getTemperature() << " K" << std::endl;
        std::cerr << "    delta t     : "
            << std::setw(7) << std::fixed << std::setprecision(3)
            << integrator_.getStepSize() / Constant::cafetime << " cafetime"
            << std::endl;

        observer_.initialize(total_step);
    }

    void initialize(const std::vector<OpenMM::Vec3>& initial_position)
    {
        // Set starting positions of the atoms.
        context_.setPositions(initial_position);
    }

    void run()
    {
        const std::size_t total_frame = std::floor(total_step_/save_step_);
        for(std::size_t frame_num=0; frame_num<total_frame; ++frame_num)
        {
            observer_.output(frame_num*save_step_, context_);

            integrator_.step(save_step_);
        }

        observer_.finalize();
    }

    const std::size_t total_step() const noexcept { return total_step_; }
    const std::size_t save_step()  const noexcept { return save_step_; }

  private:
    std::unique_ptr<OpenMM::System> system_ptr_;
    OpenMM::LangevinIntegrator      integrator_;
    OpenMM::Context                 context_;
    std::size_t                     total_step_;
    std::size_t                     save_step_;
    Observer                        observer_;
};

#endif // OPEN_AICG2_PLUS_SIMULATOR_HPP
