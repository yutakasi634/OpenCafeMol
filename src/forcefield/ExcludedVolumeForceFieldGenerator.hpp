#ifndef OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <set>
#include <optional>

class ExcludedVolumeForceFieldGenerator final: public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    // The size of the vector representing the per-particle parameters
    // (in this case, radiuses) must match the system size. This must be guaranteed
    // inside the read function.
    ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept override { return "ExcludedVolume"; }

  private:
    double                              eps_;
    double                              cutoff_;
    std::vector<std::optional<double>>  radiuses_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
