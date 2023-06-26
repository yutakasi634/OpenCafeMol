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
        const bool use_periodic, const std::size_t ffgen_id = 0)
        : indices_vec_(indices_vec), ks_(ks),
          fourier_table_(fourier_table), aa_pair_name_(aa_pair_name),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_id))
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
        const std::string fld_expression =
            "FLD"+ffgen_id_str_+"_k *"
            "(FLD"+ffgen_id_str_+"_c + FLD"+ffgen_id_str_+"_kcos1*cos(  theta) +"
            "                          FLD"+ffgen_id_str_+"_ksin1*sin(  theta) +"
            "                          FLD"+ffgen_id_str_+"_kcos2*cos(2*theta) +"
            "                          FLD"+ffgen_id_str_+"_ksin2*sin(2*theta) +"
            "                          FLD"+ffgen_id_str_+"_kcos3*cos(3*theta) +"
            "                          FLD"+ffgen_id_str_+"_ksin3*sin(3*theta))";
        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(fld_expression);
        torsion_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_k");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_c");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_kcos1");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_ksin1");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_kcos2");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_ksin2");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_kcos3");
        torsion_ff->addPerTorsionParameter("FLD"+ffgen_id_str_+"_ksin3");

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
    const std::string name() const noexcept
    {
        return "FlexibleLocalDihedral (" + aa_pair_name_ + ")";
    }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::array<double, 7>     fourier_table_;
    std::string               aa_pair_name_;
    const bool                use_periodic_;
    const std::string         ffgen_id_str_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
