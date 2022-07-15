#ifndef OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <limits>
#include <set>
#include <optional>
#include <OpenMM.h>

class ExcludedVolumeForceFieldGenerator
{
  public:
    using exclusion_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const exclusion_pairs_type& exclusion_pairs)
        : eps_(eps), cutoff_(cutoff),
          radiuses_(radiuses), exclusion_pairs_(exclusion_pairs)
    {}

    std::unique_ptr<OpenMM::CustomNonbondedForce> generate() const noexcept
    {
        std::cerr << "    Global        : ExcludedVolume" << std::endl;

        const std::string potential_formula =
            "epsilon*(sigma1+sigma2)^12*((1/r)^12-cutoff_correction)";
        auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);


        exv_ff->addPerParticleParameter("sigma");
        exv_ff->addGlobalParameter("epsilon", eps_);

        double max_radius        = std::numeric_limits<double>::min();
        double second_max_radius = std::numeric_limits<double>::min();
        const auto& itr = std::find(radiuses_.begin(), radiuses_.end(), std::nullopt);
        if(itr != radiuses_.end())
        {
            std::set<int> participants;
            for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
            {
                const std::optional<double>& radius = radiuses_[idx];
                if(radius)
                {
                    exv_ff->addParticle({radius.value()});
                    participants.insert(idx);

                    if(max_radius < radius.value())
                    {
                        max_radius        = radius.value();
                        second_max_radius = max_radius;
                    }
                }
                else
                {
                    exv_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
                }
            }
            exv_ff->addInteractionGroup(participants, participants);
        }
        else
        {
            for(const auto& radius : radiuses_)
            {
                exv_ff->addParticle({radius.value()});
                if(max_radius < radius.value())
                 {
                    max_radius        = radius.value();
                    second_max_radius = max_radius;
                 }
            }
        }

        // set cutoff
        exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        const double cutoff_distance   = (max_radius + second_max_radius)*cutoff_;
        const double cutoff_correction = std::pow(1.0 / cutoff_distance, 12);
        exv_ff->setCutoffDistance(cutoff_distance);
        exv_ff->addGlobalParameter("cutoff_correction", cutoff_correction);

        for(const auto& exclusion_pair : exclusion_pairs_)
        {
            exv_ff->addExclusion(exclusion_pair.first, exclusion_pair.second);
        }

        return exv_ff;
    }

  private:
    double                             eps_;
    double                             cutoff_;
    std::vector<std::optional<double>> radiuses_;
    exclusion_pairs_type               exclusion_pairs_;
};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
