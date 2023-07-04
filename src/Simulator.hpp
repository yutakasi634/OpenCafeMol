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
        const SystemGenerator& system_gen, const IntegratorGeneratorBase& integrator_gen,
        const std::vector<OpenMM::Vec3>& initial_position,
        const std::size_t total_step, const std::size_t save_step,
        std::vector<std::unique_ptr<ObserverBase>>& observers, bool output_progress = true)
        : Simulator(system_gen, integrator_gen,
                    total_step, save_step, observers, output_progress)
    {
        this->initialize(initial_position);
    }

    Simulator(
        const SystemGenerator& system_gen, const IntegratorGeneratorBase& integrator_gen,
        const std::size_t total_step, const std::size_t save_step,
        std::vector<std::unique_ptr<ObserverBase>>& observers, bool output_progress = true)
        : system_ptr_(system_gen.generate().release()), integrator_ptr_(integrator_gen.generate()),
          context_(*system_ptr_, *integrator_ptr_, OpenMM::Platform::getPlatformByName("CUDA")),
          total_step_(total_step), save_step_(save_step), output_progress_(output_progress),
          progress_bar_(/* width of bar = */ 50)
    {
        std::cerr << "initializing simulator..." << std::endl;
        std::cerr << "    total step : " << total_step_ << " step" << std::endl;
        std::cerr << "    save step  : " << save_step_  << " step" << std::endl;
        std::cerr << "integrator information" << std::endl;
        std::cerr << integrator_gen.dump_info();

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
            if(output_progress_)
            {
                progress_bar_.finalize(std::cerr);
            }
        }
    }

    std::size_t total_step() const noexcept { return total_step_; }
    std::size_t save_step()  const noexcept { return save_step_; }

  private:
    std::unique_ptr<OpenMM::System>            system_ptr_;
    std::unique_ptr<OpenMM::Integrator>        integrator_ptr_;
    OpenMM::Context                            context_;
    std::size_t                                total_step_;
    std::size_t                                save_step_;
    std::vector<std::unique_ptr<ObserverBase>> observers_;
    bool                                       output_progress_;
    ProgressBar                                progress_bar_;
};

#endif // OPEN_AICG2_PLUS_SIMULATOR_HPP
