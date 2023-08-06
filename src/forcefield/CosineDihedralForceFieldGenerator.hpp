#ifndef OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <regex>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

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

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        std::string potential_formula =  "k_periodic * (1 - cos(n0 * (theta - t0)))";
        // The 3SPN2C DNA model paper (Freeman et al., JCP, 2014) shows k*(1 + cos x) formula.
        // However, we should note that this is a misprint of k*(1 - cos x).
        // Indeed, the other implementations of 3SPN2 such as open3spn2, LAMMPS, and CafeMol
        // use k* (1 - cos x) formula.

        const std::map<std::string, std::string> ff_params =
        {
            {"k_periodic", ffgen_id_ + "_k_periodic"},
            {"t0",         ffgen_id_ + "_t0"},
            {"n0",         ffgen_id_ + "_n0"},
        };

        for(const auto& param : ff_params)
        {
              potential_formula = std::regex_replace(
                potential_formula, std::regex(param.first), param.second);
        }

        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(potential_formula);

        torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);

        torsion_ff->addPerTorsionParameter(ff_params.at("k_periodic"));
        torsion_ff->addPerTorsionParameter(ff_params.at("t0"));
        torsion_ff->addPerTorsionParameter(ff_params.at("n0"));

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& indices = indices_vec_[idx];
            torsion_ff->addTorsion(
                    indices[0], indices[1], indices[2], indices[3],
                    {ks_[idx], theta0s_[idx], ns_[idx]});
        }

        return torsion_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string                      name()    const noexcept { return "CosineDihedral"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       theta0s_;
    std::vector<double>       ns_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
