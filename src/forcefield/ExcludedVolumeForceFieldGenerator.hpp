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
    ExcludedVolumeForceFieldGenerator(const double eps,
        const std::vector<std::optional<double>>& radiuses,
        const exclusion_pairs_type& exclusion_pairs)
        : eps_(eps), radiuses_(radiuses), exclusion_pairs_(exclusion_pairs)
    {}

    std::unique_ptr<OpenMM::CustomNonbondedForce> generate() const noexcept
    {
        std::cerr << "    Global        : ExcludedVolume" << std::endl;

        // TODO: add cutoff
        const std::string potential_formula = "epsilon*((sigma1+sigma2)/r)^12";
        auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);
        exv_ff->addPerParticleParameter("sigma");
        exv_ff->addGlobalParameter("epsilon", eps_);

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
            }
        }

        for(const auto& exclusion_pair : exclusion_pairs_)
        {
            exv_ff->addExclusion(exclusion_pair.first, exclusion_pair.second);
        }

        return exv_ff;
    }

  private:
    double                             eps_;
    std::vector<std::optional<double>> radiuses_;
    exclusion_pairs_type               exclusion_pairs_;
};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
