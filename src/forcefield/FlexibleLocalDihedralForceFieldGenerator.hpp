#ifndef OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

class FlexibleLocalDihedralForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::array<std::size_t, 4>;

  public:
    FlexibleLocalDihedralForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::array<double, 7>& fourier_table, const std::string& aa_pair_name,
        const bool use_periodic)
        : indices_vec_(indices_vec), ks_(ks),
          fourier_table_(fourier_table), aa_pair_name_(aa_pair_name),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("FLD{}", ffid.gen()))
    {
        if(!(indices_vec.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] FlexibleLocalDihedralForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

    std::string name() const noexcept override
    {
        return "FlexibleLocalDihedral (" + aa_pair_name_ + ")";
    }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::array<double, 7>     fourier_table_;
    std::string               aa_pair_name_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
