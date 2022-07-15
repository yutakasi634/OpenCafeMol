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
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& bonded_pairs,
        const index_pairs_type& additional_exclusion_pairs)
        : eps_(eps), cutoff_(cutoff), radiuses_(radiuses),
          bonded_pairs_(bonded_pairs),
          additional_exclusion_pairs_(additional_exclusion_pairs)
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

        create_exclusion_from_bond(exv_ff, bonded_pairs_);
        add_exclusion(exv_ff, additional_exclusion_pairs_);

        return exv_ff;
    }

  private:
    void create_exclusion_from_bond(
            std::unique_ptr<OpenMM::CustomNonbondedForce>& exv_ff,
            const index_pairs_type& bonded_pairs) const noexcept
    {
        exv_ff->createExclusionsFromBonds({bonded_pairs.begin(), bonded_pairs.end()}, 3);
    }

    void add_exclusion(
            std::unique_ptr<OpenMM::CustomNonbondedForce>& exv_ff,
            const index_pairs_type& additional_exclusion_pairs) const noexcept
    {
        const std::size_t num_exclusions = exv_ff->getNumExclusions();

        std::vector<std::pair<std::size_t, std::size_t>> exclusion_pairs;
        for(std::size_t pair_idx=0; pair_idx<num_exclusions; ++pair_idx)
        {
            int first, second;
            exv_ff->getExclusionParticles(pair_idx, first, second);
            exclusion_pairs.push_back(std::make_pair(first, second));
        }

        for(const auto& pair : additional_exclusion_pairs)
        {
            const auto result =
                std::find(exclusion_pairs.begin(), exclusion_pairs.end(), pair);
            if(result == exclusion_pairs.end())
            {
                exv_ff->addExclusion(pair.first, pair.second);
            }
        }

    }

  private:
    double                             eps_;
    double                             cutoff_;
    std::vector<std::optional<double>> radiuses_;
    const index_pairs_type&            bonded_pairs_;
    const index_pairs_type&            additional_exclusion_pairs_;

};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
