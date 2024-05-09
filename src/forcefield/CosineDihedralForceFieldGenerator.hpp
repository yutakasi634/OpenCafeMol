#ifndef OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class CosineDihedralForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type    = std::array<std::size_t, 4>;
    using index_pair_type = std::pair<std::size_t, std::size_t>;

  public:
    CosineDihedralForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::vector<double>&       theta0s,     const std::vector<double>& ns,
        const bool use_periodic)
        : indices_vec_(indices_vec), ks_(ks), theta0s_(theta0s), ns_(ns),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("CD{}", ffid.gen()))
    {
        const std::size_t system_size = indices_vec.size();
        if(!(system_size == ks.size() && system_size == theta0s.size() &&
             system_size == ns.size()))
        {
            std::ostringstream oss;
            oss << "[error] CosineDihedralForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << ks.size()          << "), "
                   "theta0 ("      << theta0s.size()     << ") and "
                   "n ("           << ns.size()          << " is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string                      name()    const override { return "CosineDihedral"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       theta0s_;
    std::vector<double>       ns_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
