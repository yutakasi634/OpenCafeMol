#ifndef OPEN_AICG2_PLUS_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>

class HarmonicBondForceFieldGenerator
{
  public:
    using index_pair_type = std::pair<std::size_t, std::size_t>;

  public:
    HarmonicBondForceFieldGenerator(
        const std::vector<index_pair_type> indices_vec,
        const std::vector<double>& v0s, const std::vector<double>& ks)
        : indices_vec_(indices_vec), v0s_(v0s), ks_(ks)
    {
        if(!(indices_vec.size() == v0s.size() && v0s.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] HarmonicBondForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "v0 ("          << v0s.size()         << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                << "The number os these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::HarmonicBondForce> generate() const noexcept
    {
        auto bond_ff = std::make_unique<OpenMM::HarmonicBondForce>();
        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const index_pair_type& idx_pair = indices_vec_[idx];
            bond_ff->addBond(idx_pair.first, idx_pair.second, v0s_[idx], ks_[idx]);
        }

        return bond_ff;
    }

    void add_exclusion(
            std::vector<index_pair_type>& exclusion_pairs) const noexcept
    {
        for(const auto& indices : indices_vec_)
        {
            exclusion_pairs.push_back(std::make_pair(indices.first, indices.second));
        }
    }

  private:
    std::vector<index_pair_type> indices_vec_;
    std::vector<double>          v0s_;
    std::vector<double>          ks_;
};

#endif // OPEN_AICG2_PLUS_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP
