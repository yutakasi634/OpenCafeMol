#ifndef OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

class UniformWeeksChandlerAndersenForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    UniformWeeksChandlerAndersenForceFieldGenerator(
        const std::size_t system_size, const double eps, const double sigma,
        const std::vector<std::size_t> former_group_vec,
        const std::vector<std::size_t> latter_group_vec,
        const index_pairs_type& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::size_t former_group_size() const noexcept { return former_group_size_; }
    std::size_t latter_group_size() const noexcept { return latter_group_size_; }

    std::string name() const noexcept override { return "UniformWeeksChadlerAndersen"; }

  private:
    std::size_t      system_size_;
    double           eps_;
    double           sigma_;
    index_pairs_type ignore_list_;
    bool             use_periodic_;
    std::string      ffgen_id_;
    std::size_t      former_group_size_;
    std::size_t      latter_group_size_;

    std::vector<interaction_group_type> interaction_groups_;
};

#endif // OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
