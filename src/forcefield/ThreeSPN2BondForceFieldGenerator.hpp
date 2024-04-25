#ifndef OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

class ThreeSPN2BondForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    ThreeSPN2BondForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& k2s, const std::vector<double>& k4s, const std::vector<double>& v0s,
        const bool use_periodic)
        : indices_vec_(indices_vec), k2s_(k2s), k4s_(k4s), v0s_(v0s),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("TSPN2B{}", ffid.gen()))
    {
        if(!(indices_vec.size() == v0s.size() && v0s.size() == k2s.size()))
        {
            std::ostringstream oss;
            oss << "[error] ThreeSPN2BondForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << k2s.size()         << "), "
                   "v0 ("          << v0s.size()         << ") is not matched "
               << "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        std::string potential_formula = fmt::format(
            "{id}_k2 * (r - {id}_v0)^2 + 100.0 * {id}_k4 * (r - {id}_v0)^4",
            fmt::arg("id", ffgen_id_));

        auto bond_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);

        bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        bond_ff->addPerBondParameter(fmt::format("{}_k2", ffgen_id_));
        bond_ff->addPerBondParameter(fmt::format("{}_k4", ffgen_id_));
        bond_ff->addPerBondParameter(fmt::format("{}_v0", ffgen_id_));

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& idx_pair = indices_vec_[idx];
            bond_ff->addBond(idx_pair.first, idx_pair.second,
                             {k2s_[idx], k4s_[idx], v0s_[idx]});
        }

        return bond_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string                      name()    const noexcept final { return "3SPN2Bond"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       k2s_;
    std::vector<double>       k4s_;
    std::vector<double>       v0s_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP
