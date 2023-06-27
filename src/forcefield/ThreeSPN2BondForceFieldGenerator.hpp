#ifndef OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP

#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>

class ThreeSPN2BondForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    ThreeSPN2BondForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& k2s, const std::vector<double>& k4s, const std::vector<double>& v0s,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : indices_vec_(indices_vec), k2s_(k2s), k4s_(k4s), v0s_(v0s),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
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
        std::string potential_formula = "k2 * (r - v0)^2 + 100.0 * k4 * (r - v0)^4";

        const std::map<std::string, std::string> ff_params =
        {
            {"k2", "TSPN2B" + ffgen_id_str_ + "_k2"},
            {"k4", "TSPN2B" + ffgen_id_str_ + "_k4"},
            {"v0", "TSPN2B" + ffgen_id_str_ + "_v0"},
        };
        for(auto itr = ff_params.begin(); itr != ff_params.end(); ++itr)
        {
              potential_formula = std::regex_replace(
                potential_formula, std::regex(itr->first), itr->second);
        }

        auto bond_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);

        bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        bond_ff->addPerBondParameter(ff_params.at("k2"));
        bond_ff->addPerBondParameter(ff_params.at("k4"));
        bond_ff->addPerBondParameter(ff_params.at("v0"));

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& idx_pair = indices_vec_[idx];
            bond_ff->addBond(idx_pair.first, idx_pair.second,
                             {k2s_[idx], k4s_[idx], v0s_[idx]});
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

    const std::vector<indices_type>& indices() const noexcept {return indices_vec_; }
    std::string name() const noexcept {return "3SPN2Bond";};

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       k2s_;
    std::vector<double>       k4s_;
    std::vector<double>       v0s_;
    bool                      use_periodic_;
    std::string               ffgen_id_str_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_BOND_FORCE_FIELD_GENERATOR_HPP
