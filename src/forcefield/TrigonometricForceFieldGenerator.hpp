#ifndef OPEN_AICG2_PLUS_TRIGONOMETRIC_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_TRIGONOMETRIC_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <optional>

class TrigonometricForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    TrigonometricForceFieldGenerator(
        const std::vector<std::optional<double>>& epsilons,
        const std::vector<std::optional<double>>& sigmas,
        const std::vector<std::optional<double>>& omegas,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override { return "Trigonometric"; }

  private:
    std::vector<std::optional<double>>  epsilons_;
    std::vector<std::optional<double>>  sigmas_;
    std::vector<std::optional<double>>  omegas_;
    index_pairs_type                    ignore_list_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;
    std::vector<interaction_group_type> interaction_groups_;
};

#endif // OPEN_AICG2_PLUS_TRIGONOMETRIC_FORCE_FIELD_GENERATOR_HPP
