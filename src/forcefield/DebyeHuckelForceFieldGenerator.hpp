#ifndef OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_DEBYE_HUCKEL_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <optional>
#include <set>
#include <utility>
#include <vector>

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
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override { return "DebyeHuckel"; }

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
