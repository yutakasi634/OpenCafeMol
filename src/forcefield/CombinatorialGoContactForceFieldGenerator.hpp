#ifndef OPEN_AICG2_PLUS_COMBINATORIAL_GO_CONTACT_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_COMBINATORIAL_GO_CONTACT_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

class CombinatorialGoContactForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    CombinatorialGoContactForceFieldGenerator(
        const double k, const double r0, const double cutoff_ratio,
        const std::vector<std::size_t> donor_indices_vec,
        const std::vector<std::size_t> acceptor_indices_vec,
        const index_pairs_type& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : k_(k), r0_(r0), cutoff_ratio_(cutoff_ratio),
          donor_indices_vec_(donor_indices_vec),
          acceptor_indices_vec_(acceptor_indices_vec),
          ignore_list_(ignore_list), use_periodic_(use_periodic),
          ffgen_id_(fmt::format("CGC{}", ffid.gen()))
    {}

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override { return "CombinatorialGoContact"; }

  private:
    double                   k_;
    double                   r0_;
    double                   cutoff_ratio_;
    std::vector<std::size_t> donor_indices_vec_;
    std::vector<std::size_t> acceptor_indices_vec_;
    index_pairs_type         ignore_list_;
    bool                     use_periodic_;
    std::string              ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_COMBINATORIAL_GO_CONTACT_FORCE_FIELD_GENERATOR_HPP
