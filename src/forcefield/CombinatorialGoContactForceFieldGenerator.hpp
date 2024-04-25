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

    std::unique_ptr<OpenMM::Force> generate() const override
    {
        std::string potential_formula = fmt::format(
            "{id}_k *"
                "(5 * ({id}_r0 / distance(a1, d1))^12 -"
                "6 * ({id}_r0 / distance(a1, d1))^10)",
            fmt::arg("id", ffgen_id_));

        auto contact_ff = std::make_unique<OpenMM::CustomHbondForce>(potential_formula);

        contact_ff->addGlobalParameter(fmt::format("{}_k",  ffgen_id_), k_);
        contact_ff->addGlobalParameter(fmt::format("{}_r0", ffgen_id_), r0_);

        // set cutoff
        if(use_periodic_)
        {
            contact_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
        }
        else
        {
            contact_ff->setNonbondedMethod(OpenMM::CustomHbondForce::CutoffNonPeriodic);
        }
        contact_ff->setCutoffDistance(r0_ * cutoff_ratio_);

        for(const auto& donor_particles : donor_indices_vec_)
        {
            contact_ff->addDonor(donor_particles, -1, -1);
        }

        for(const auto& acceptor_particles : acceptor_indices_vec_)
        {
            contact_ff->addAcceptor(acceptor_particles, -1, -1);
        }

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            contact_ff->addExclusion(pair.first, pair.second);
        }

        return contact_ff;
    }

    std::string name() const noexcept { return "CombinatorialGoContact"; }

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
