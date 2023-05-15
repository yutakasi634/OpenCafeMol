#ifndef OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"

class GoContactForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    GoContactForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& ks, const std::vector<double>& r0s,
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : indices_vec_(indices_vec), ks_(ks), r0s_(r0s),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
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

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula =
            "GC"+ffgen_id_str_+"_k*"
            "(5*(GC"+ffgen_id_str_+"_r0/r)^12-6*(GC"+ffgen_id_str_+"_r0/r)^10)";
        auto contact_ff = std::make_unique<OpenMM::CustomBondForce>(potential_formula);
        contact_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        contact_ff->addPerBondParameter("GC"+ffgen_id_str_+"_k");
        contact_ff->addPerBondParameter("GC"+ffgen_id_str_+"_r0");

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const indices_type& indices = indices_vec_[idx];
            contact_ff->addBond(indices.first, indices.second, {ks_[idx], r0s_[idx]});
        }

        return contact_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    const std::string name() const noexcept { return "GoContact"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       r0s_;
    const bool                use_periodic_;
    const std::string         ffgen_id_str_;
};

#endif // OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
