#ifndef OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include "src/util/Constants.hpp"

class DebyeHuckelForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    DebyeHuckelForceFieldGenerator(const double ionic_strength,
        const double temperature, const double cutoff_ratio,
        const std::vector<std::optional<double>>& charges,
        const index_pairs_type& ignore_list)
        : ionic_strength_(ionic_strength), temperature_(temperature),
          cutoff_ratio_(cutoff_ratio), charges_(charges), ignore_list_(ignore_list)
    {
        const double epsk    = calc_dielectric_water(temperature_, ionic_strength_);
        const double eps0_ee = Constant::eps0 / Constant::elementary_charge
                                              / Constant::elementary_charge;
        inv_4_pi_eps0_epsk_ = 1.0 / (4 * Constant::pi * eps0_ee * epsk);

        // convert [M] (mol/L) to [mol/m^3]
        const double I = ionic_strength_ * 1.0e-3; // [M/m^3]

        debye_length_ = std::sqrt((eps0_ee * epsk * Constant::kB * temperature_) /
                                 (2. * Constant::Na * I)) * 10e9; // [nm]
        abs_cutoff_   = debye_length_ * cutoff_ratio_;
        cutoff_correction_ = std::exp(cutoff_ratio_) / abs_cutoff_;
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {

        const std::string potential_formula =
            "q1*q2*inv_4_pi_eps0_epsk*(exp(-r/debye_length)-cutoff_correction)";
        auto dh_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        dh_ff->addPerParticleParameter("q");
        dh_ff->addGlobalParameter("inv_4_pi_eps0_epsk", inv_4_pi_eps0_epsk_);
        dh_ff->addGlobalParameter("debye_length",       debye_length_);
        dh_ff->addGlobalParameter("cutoff_correction",  cutoff_correction_);

        const auto& itr = std::find(charges_.begin(), charges_.end(), std::nullopt);
        if(itr != charges_.end())
        {
            std::set<int> participants;
            for(std::size_t idx=0; idx<charges_.size(); ++idx)
            {
                const std::optional<double>& charge = charges_[idx];
                if(charge)
                {
                    dh_ff->addParticle({charge.value()});
                    participants.insert(idx);
                }
                else
                {
                    dh_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
                }
            }
            dh_ff->addInteractionGroup(participants, participants);
        }
        else
        {
            // all system particles are participants case
            for(const auto& charge : charges_)
            {
                dh_ff->addParticle({charge.value()});
            }
        }

        // set cutoff
        dh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        dh_ff->setCutoffDistance(abs_cutoff_);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            dh_ff->addExclusion(pair.first, pair.second);
        }

        return dh_ff;
    }

  private:
    double calc_dielectric_water(const double T, const double C) const noexcept
    {
        return (249.4 - 0.778 * T + 7.2e-4 * T * T) *
            (1. - 2.551e-1 * C + 5.151e-2 * C * C - 6.889e-3 * C * C * C);
    }

  private:
    const double                             ionic_strength_; // [M]
    const double                             temperature_;    // [K]
    const double                             cutoff_ratio_;   // relative to the debye length
    const std::vector<std::optional<double>> charges_;
    index_pairs_type                         ignore_list_;

    double debye_length_;
    double inv_4_pi_eps0_epsk_;
    double cutoff_correction_;
    double abs_cutoff_;
};

#endif // OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP
