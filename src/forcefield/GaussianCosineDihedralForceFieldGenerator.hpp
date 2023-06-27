#ifndef OPEN_AICG2_PLUS_GAUSSIAN_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GAUSSIAN_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include <regex>
#include <iostream> // debug

#include "src/util/Constants.hpp"

class GaussianCosineDihedralForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type    = std::array<std::size_t, 4>;
    using index_pair_type = std::pair<std::size_t, std::size_t>;

  public:
    GaussianCosineDihedralForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& ks_gauss,      const std::vector<double>& ks_cos,
        const std::vector<double>& theta0s_gauss, const std::vector<double>& theta0s_cos,
        const std::vector<double>& sigmas,        const std::vector<double>& ns,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : indices_vec_(indices_vec), ks_gauss_(ks_gauss), ks_cos_(ks_cos),
          theta0s_gauss_(theta0s_gauss), theta0s_cos_(theta0s_cos),
          sigmas_(sigmas), ns_(ns),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
    {
        const std::size_t system_size = indices_vec.size();
        if(!(system_size == ks_gauss.size() && system_size == ks_cos.size() &&
             system_size == theta0s_gauss.size() && system_size == theta0s_cos.size() &&
             system_size == sigmas.size() && system_size == ns.size()))
        {
            std::ostringstream oss;
            oss << "[error] GaussianCosineDihedralForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec ("  << indices_vec.size()   << "), "
                   "k_gauss ("      << ks_gauss.size()      << "), "
                   "k_cos ("        << ks_cos.size()        << "), "
                   "theta0_gauss (" << theta0s_gauss.size() << "), "
                   "theta0_cos ("   << theta0s_cos.size()   << "), "
                   "sigma ("        << sigmas.size()        << ") and, "
                   "n ("            << ns.size()            << " is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        std::string potential_formula = "energy;"
            "energy      = k_gaussian*exp(-dt_periodic^2/(2*sigma^2)) + k_cos*(1 - cs);"
            "cs          = cos(n0 * dt_cos);"
            "dt_periodic = dt_gaussian - floor((dt_gaussian + pi)/(2*pi))*(2*pi);"
            "dt_gaussian = theta - t0_gaussian;"
            "dt_cos      = theta - t0_cos;";

        // The 3SPN2C DNA model paper (Freeman et al., JCP, 2014) shows cosine term as k*(1 + cos x).
        // However, we should note that this is a misprint of k*(1 - cos x).
        // Indeed, the other implementations of 3SPN2 such as open3spn2, LAMMPS, and CafeMol
        // use cosine term as k* (1 - cos x) formula.

        const std::map<std::string, std::string> ff_params =
        {
            {"k_gaussian",  "GCD" + ffgen_id_str_ + "_k_gaussian"},
            {"k_cos",       "GCD" + ffgen_id_str_ + "_k_periodic"},
            {"t0_gaussian", "GCD" + ffgen_id_str_ + "_t0_gauss"},
            {"t0_cos",      "GCD" + ffgen_id_str_ + "_t0_cos"},
            {"sigma",       "GCD" + ffgen_id_str_ + "_sigma"},
            {"n0",          "GCD" + ffgen_id_str_ + "_n0"},
        };

        for(auto itr = ff_params.begin(); itr != ff_params.end(); ++itr)
        {
              potential_formula = std::regex_replace(
                potential_formula, std::regex(itr->first), itr->second);
        }

        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(potential_formula);

        torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);

        torsion_ff->addPerTorsionParameter(ff_params.at("k_gaussian"));
        torsion_ff->addPerTorsionParameter(ff_params.at("k_cos"));
        torsion_ff->addPerTorsionParameter(ff_params.at("t0_gaussian"));
        torsion_ff->addPerTorsionParameter(ff_params.at("t0_cos"));
        torsion_ff->addPerTorsionParameter(ff_params.at("sigma"));
        torsion_ff->addPerTorsionParameter(ff_params.at("n0"));
        torsion_ff->addGlobalParameter    ("pi", Constant::pi);

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& indices = indices_vec_[idx];
            torsion_ff->addTorsion(
                    indices[0], indices[1], indices[2], indices[3],
                    {ks_gauss_[idx], ks_cos_[idx], theta0s_gauss_[idx],
                     theta0s_cos_[idx], sigmas_[idx], ns_[idx]});
        }

        return torsion_ff;
    }

    void add_exclusion(std::vector<index_pair_type>& exclusion_pairs) const noexcept
    {
        for(const auto& indices : indices_vec_)
        {
            exclusion_pairs.push_back(std::make_pair(indices[0], indices[3]));
        }
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string name() const noexcept { return "PeriodicGaussianCosineDihedral"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_gauss_;
    std::vector<double>       ks_cos_;
    std::vector<double>       theta0s_gauss_;
    std::vector<double>       theta0s_cos_;
    std::vector<double>       sigmas_;
    std::vector<double>       ns_;
    const bool                use_periodic_;
    const std::string         ffgen_id_str_;
};

#endif // OPEN_AICG2_PLUS_GAUSSIAN_COSINE_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
