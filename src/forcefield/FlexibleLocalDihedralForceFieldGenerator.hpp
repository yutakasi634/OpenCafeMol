#ifndef OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"

class FlexibleLocalDihedralForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::array<std::size_t, 4>;

  public:
    FlexibleLocalDihedralForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::array<double, 7>& fourier_table, const std::string& aa_pair_name,
        const bool use_periodic, const std::size_t ffgen_id)
        : indices_vec_(indices_vec), ks_(ks),
          fourier_table_(fourier_table), aa_pair_name_(aa_pair_name),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("FLD{}", ffgen_id))
    {
        if(!(indices_vec.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] FlexibleLocalDihedralForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string fld_expression = fmt::format("{id}_k * ("
                "{id}_c +"
                "{id}_kcos1*cos(  theta) +"
                "{id}_ksin1*sin(  theta) +"
                "{id}_kcos2*cos(2*theta) +"
                "{id}_ksin2*sin(2*theta) +"
                "{id}_kcos3*cos(3*theta) +"
                "{id}_ksin3*sin(3*theta)"
            ")", fmt::arg("id", ffgen_id_));

        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(fld_expression);
        torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        torsion_ff->addPerTorsionParameter(fmt::format("{}_k",     ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_c",     ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_kcos1", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_ksin1", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_kcos2", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_ksin2", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_kcos3", ffgen_id_));
        torsion_ff->addPerTorsionParameter(fmt::format("{}_ksin3", ffgen_id_));

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const auto&  indices = indices_vec_[idx];
            const double k       = ks_[idx];
            torsion_ff->addTorsion(
                    indices[0], indices[1], indices[2], indices[3],
                    {k, fourier_table_[0],
                     fourier_table_[1], fourier_table_[2], fourier_table_[3],
                     fourier_table_[4], fourier_table_[5], fourier_table_[6]});
        }

        return torsion_ff;
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string name() const noexcept
    {
        return "FlexibleLocalDihedral (" + aa_pair_name_ + ")";
    }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::array<double, 7>     fourier_table_;
    std::string               aa_pair_name_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
