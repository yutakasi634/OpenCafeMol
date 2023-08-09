#ifndef OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>

#include <fmt/core.h>

#include "src/util/Constants.hpp"
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

class DebyeHuckelForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    DebyeHuckelForceFieldGenerator(const double ionic_strength,
        const double temperature, const double cutoff_ratio,
        const std::vector<std::optional<double>>& charges,
        const index_pairs_type& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : ionic_strength_(ionic_strength), temperature_(temperature),
          cutoff_ratio_(cutoff_ratio), charges_(charges), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("DH{}", ffid.gen()))
    {
        const double epsk = calc_dielectric_water(temperature_, ionic_strength_); // dimensionless
        const double eps0 = Constant::eps0 / Constant::elementary_charge
                                           / Constant::elementary_charge; // [mol/KJ/nm]
        std::cerr << "        eps0 is " << eps0 << " mol/KJ/nm" << std::endl;
        std::cerr << "        epsk is " << epsk << " (dimensionless)" << std::endl;

        inv_4_pi_eps0_epsk_ = 1.0 / (4.0 * Constant::pi * eps0 * epsk); // [KJ nm/mol]

        // convert [M] (mol/L) to [mol/nm^3]
        const double I = ionic_strength_ * 1.0e-24/*[/L]->[/nm^3]*/; // [mol/nm^3]

        debye_length_ =
            std::sqrt((eps0 * epsk * Constant::kB * 1.0e-3 /*[J]->[KJ]*/ * temperature_) / (2. * I)); // [nm]
        std::cerr << "        debye length is " << debye_length_ << " nm" << std::endl;
        abs_cutoff_   = debye_length_ * cutoff_ratio_;
        cutoff_correction_ = std::exp(-cutoff_ratio_) / abs_cutoff_;

        if(ignore_group_pairs.size() == 0)
        {
            // without interaction group force calculation is faster than
            // with that. So if all particle have the parameter for this
            // calculation, force don't use interaction group
            if(std::any_of(charges_.begin(), charges_.end(),
                          [](const std::optional<double>& val) { return !val; }))
            {
                std::set<int> participants;
                for(std::size_t idx=0; idx<charges_.size(); ++idx)
                {
                    if(charges_[idx])
                    {
                        participants.insert(idx);
                    }
                }
                interaction_groups_.push_back({ participants, participants });
            }
            else
            {
                std::cerr << "        all particles are participants in this interaction" << std::endl;
            }
        }
        else // make interaction group when group based ignoration specified
        {
            std::set<std::string> related_group_names;
            for(const auto& name_group_pair : ignore_group_pairs)
            {
                related_group_names.insert(name_group_pair.first);
                related_group_names.insert(name_group_pair.second);
            }

            std::map<std::string, std::set<int>> related_group_map;
            for(const auto& name : related_group_names)
            {
                related_group_map.insert(std::make_pair(name, std::set<int>()));
            }

            std::set<int> others;
            for(std::size_t idx=0; idx<charges_.size(); ++idx)
            {
                if(charges_[idx])
                {
                    if(group_vec[idx])
                    {
                        const std::string group_name = group_vec[idx].value();
                        if(related_group_names.count(group_name) != 0)
                        {
                            related_group_map.at(group_name).insert(idx);
                        }
                        else
                        {
                            others.insert(idx);
                        }
                    }
                    else
                    {
                        others.insert(idx);
                    }
                }
            }

            std::vector<std::pair<std::string, std::set<int>>> related_group_vec;
            for(const auto& name_group_pair : related_group_map)
            {
                related_group_vec.push_back(name_group_pair);
            }

            interaction_groups_.push_back({ others, others });
            for(std::size_t idx_i=0 ; idx_i<related_group_vec.size(); ++idx_i)
            {
                const auto& name_group_pair_i = related_group_vec[idx_i];
                const std::string&   first_name  = name_group_pair_i.first;
                const std::set<int>& first_group = name_group_pair_i.second;
                interaction_groups_.push_back({ first_group, others });
                for(std::size_t idx_j=idx_i; idx_j<related_group_vec.size(); ++idx_j)
                {
                    const auto& name_group_pair_j = related_group_vec[idx_j];
                    const std::string&   second_name  = name_group_pair_j.first;
                    if(!Utility::contains(ignore_group_pairs, { first_name, second_name }))
                    {
                        const std::set<int>& second_group = name_group_pair_j.second;
                        interaction_groups_.push_back({ first_group, second_group });
                    }
                }
            }
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula = fmt::format(
            "{id}_q1 * {id}_q2 * {id}_inv_4_pi_eps0_epsk * ("
                "exp(-r / {id}_debye_length) / r - {id}_cutoff_correction"
            ")", fmt::arg("id", this->ffgen_id_));

        auto dh_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        dh_ff->addPerParticleParameter(fmt::format("{}_q", ffgen_id_));
        dh_ff->addGlobalParameter(fmt::format("{}_inv_4_pi_eps0_epsk", ffgen_id_), inv_4_pi_eps0_epsk_); // [KJ nm /mol]
        dh_ff->addGlobalParameter(fmt::format("{}_debye_length", ffgen_id_),       debye_length_); // [nm]
        dh_ff->addGlobalParameter(fmt::format("{}_cutoff_correction", ffgen_id_),  cutoff_correction_);

        for(std::size_t idx=0; idx<charges_.size(); ++idx)
        {
            const std::optional<double>& charge = charges_[idx];
            if(charge)
            {
                dh_ff->addParticle({charge.value()});
            }
            else
            {
                dh_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
            }
        }

        // if interaction_groups size is 0, no interaction group will be added,
        // so all the particle inthe system will be considerd as participant
        for(const auto& group_pair : interaction_groups_)
        {
            dh_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            dh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            dh_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        std::cerr << "        cutoff disntace is " << abs_cutoff_ << " nm" << std::endl;
        dh_ff->setCutoffDistance(abs_cutoff_);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            dh_ff->addExclusion(pair.first, pair.second);
        }

        return dh_ff;
    }

    std::string name() const noexcept { return "DebyeHuckel"; }

  private:
    double calc_dielectric_water(const double T, const double C) const noexcept
    {
        // TODO: check this formula
        return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
            (1. - 2.551e-1 * C + 5.151e-2 * C * C - 6.889e-3 * C * C * C);
    }

  private:
    double                              ionic_strength_; // [M]
    double                              temperature_;    // [K]
    double                              cutoff_ratio_;   // relative to the debye length
    std::vector<std::optional<double>>  charges_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;

    double debye_length_;
    double inv_4_pi_eps0_epsk_;
    double cutoff_correction_;
    double abs_cutoff_;
};

#endif // OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP
