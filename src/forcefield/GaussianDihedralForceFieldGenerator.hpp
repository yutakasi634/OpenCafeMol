#ifndef OPEN_AICG2_PLUS_GAUSSIAN_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GAUSSIAN_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include <fmt/core.h>
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

class GaussianDihedralForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type    = std::array<std::size_t, 4>;
    using index_pair_type = std::pair<std::size_t, std::size_t>;

  public:
    GaussianDihedralForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::vector<double>&       theta0s,     const std::vector<double>& sigmas,
        const bool use_periodic)
        : indices_vec_(indices_vec), ks_(ks), theta0s_(theta0s), sigmas_(sigmas),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("GD{}", ffid.gen()))
    {
        const std::size_t system_size = indices_vec.size();
        if(!(system_size == ks.size() && system_size == theta0s.size() &&
             system_size == sigmas.size()))
        {
            std::ostringstream oss;
            oss << "[error] GaussianDihedralForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << ks.size()          << "), "
                   "theta0 ("      << theta0s.size()     << ") and "
                   "simga ("       << sigmas.size()      << " is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula = fmt::format(
            "{id}_k * exp(-(dt_periodic)^2 /(2*{id}_sigma^2));"
            "dt_periodic = dt - floor((dt + pi)/(2*pi))*(2*pi);"
            "dt = theta-{id}_theta0;"
            "pi = 3.1415926535897932385",
            fmt::arg("id", ffgen_id_));

        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(potential_formula);
        torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        torsion_ff->addPerTorsionParameter(fmt::format("{}_k", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_theta0", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_sigma", ffgen_id_));

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& indices = indices_vec_[idx];
            torsion_ff->addTorsion(
                    indices[0], indices[1], indices[2], indices[3],
                    {ks_[idx], theta0s_[idx], sigmas_[idx]});
        }

        return torsion_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string name() const noexcept { return "GaussianDihedral"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       theta0s_;
    std::vector<double>       sigmas_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_GAUSSIAN_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
