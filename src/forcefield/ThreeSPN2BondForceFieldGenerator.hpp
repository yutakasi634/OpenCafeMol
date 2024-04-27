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

    std::unique_ptr<OpenMM::Force> generate() const override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string                      name()    const final { return "3SPN2Bond"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       k2s_;
    std::vector<double>       k4s_;
    std::vector<double>       v0s_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP
