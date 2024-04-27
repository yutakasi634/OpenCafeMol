#ifndef OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class ThreeSPN2ExcludedVolumeForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    ThreeSPN2ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {});

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept override { return "3SPN2ExcludedVolume"; }

  private:
    double                              eps_;
    double                              cutoff_;
    std::vector<std::optional<double>>  radiuses_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;
};


// ----------------------------------------------------------------------------
// 3SPN2 parameter set
//
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
//   - TABLE III (sigma)
//   - TABLE IV  (epsilon)

struct ThreeSPN2ExcludedVolumePotentialParameter
{
    inline static const double epsilon = double(1.0); // [kJ/mol]

    inline static const std::map<std::string, double> sigma = { // [nm]
        {"P", double(4.5) * OpenMM::NmPerAngstrom},
        {"S", double(6.2) * OpenMM::NmPerAngstrom},
        {"A", double(5.4) * OpenMM::NmPerAngstrom},
        {"T", double(7.1) * OpenMM::NmPerAngstrom},
        {"G", double(4.9) * OpenMM::NmPerAngstrom},
        {"C", double(6.4) * OpenMM::NmPerAngstrom},
    };
};

#endif // OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
