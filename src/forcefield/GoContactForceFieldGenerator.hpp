#ifndef OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>

class GoContactForceFieldGenerator
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    GoContactForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& ks, const std::vector<double>& r0s)
        : indices_vec_(indices_vec), ks_(ks), r0s_(r0s)
    {
        if(!(indices_vec.size() == ks.size() && ks.size() == r0s.size()))
        {
            std::ostringstream oss;
            oss << "[error] GoContactForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << ks.size()          << ") and "
                   "r0 ("          << r0s.size()         << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::CustomBondForce> generate() const noexcept
    {
        std::cerr << "    BondLength    : GoContact" << std::endl;

        const std::string potential_formula = "k*(5*(r0/r)^12-6*(r0/r)^10)";
        auto contact_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
        contact_ff->addPerBondParameter("k");
        contact_ff->addPerBondParameter("r0");

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& indices = indices_vec_[idx];
            contact_ff->addBond(indices.first, indices.second, {ks_[idx], r0s_[idx]});
        }

        return contact_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       r0s_;
};

#endif // OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
