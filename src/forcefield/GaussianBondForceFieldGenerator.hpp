#ifndef OPEN_AICG2_PLUS_GAUSSIAN_BOND_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GAUSSIAN_BOND_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>

class GaussianBondForceFieldGenerator
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    GaussianBondForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::vector<double> v0s, const std::vector<double>& sigmas)
        : indices_vec_(indices_vec), ks_(ks), v0s_(v0s), sigmas_(sigmas)
    {
        if(!(indices_vec.size() == v0s.size() && v0s.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] GaussianBondForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << ks.size()          << "), "
                   "v0 ("          << v0s.size()         << ") and "
                   "sigmas ("      << sigmas.size()      << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::CustomBondForce> generate() const noexcept
    {
        std::cerr << "    BondLength    : Gaussian" << std::endl;

        const std::string potential_formula = "k*exp(-(r-v0)^2/(2*sigma^2))";
        auto bond_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
        bond_ff->addPerBondParameter("k");
        bond_ff->addPerBondParameter("v0");
        bond_ff->addPerBondParameter("sigma");

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& idx_pair = indices_vec_[idx];
            bond_ff->addBond(idx_pair.first, idx_pair.second,
                             {ks_[idx], v0s_[idx], sigmas_[idx]});
        }

        return bond_ff;
    }

    void add_exclusion(std::vector<indices_type>& exclusion_pairs) const noexcept
    {
        for(const auto& indices : indices_vec_)
        {
            exclusion_pairs.push_back(indices);
        }
    }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       v0s_;
    std::vector<double>       sigmas_;
};

#endif // OPEN_AICG2_PLUS_GAUSSIAN_BOND_FORCE_FIELD_GENERATOR_HPP
