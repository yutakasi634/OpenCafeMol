#ifndef OPEN_AICG2_PLUS_LENNARD_JONES_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_LENNARD_JONES_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <optional>
#include <set>
#include <utility>
#include <vector>

class LennardJonesAttractiveForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    LennardJonesAttractiveForceFieldGenerator(const double cutoff_ratio,
        const std::vector<std::optional<double>> epsilons,
        const std::vector<std::optional<double>> sigmas,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const noexcept override { return "LennardJonesAttractive"; }

  private:
    double                              cutoff_ratio_;
    std::vector<std::optional<double>>  epsilons_;
    std::vector<std::optional<double>>  sigmas_;
    index_pairs_type                    ignore_list_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;
    std::vector<interaction_group_type> interaction_groups_;
    double                              cutoff_correction_;
};

#endif // OPEN_AICG2_PLUS_LENNARD_JONES_ATTRACTIVE_FORCE_FIELD_GENERATOR_HPP
