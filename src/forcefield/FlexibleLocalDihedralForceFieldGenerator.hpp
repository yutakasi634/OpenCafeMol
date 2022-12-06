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
        const std::array<double, 7>& fourier_table, const std::string& aa_pair_name)
        : indices_vec_(indices_vec), ks_(ks),
          fourier_table_(fourier_table), aa_pair_name_(aa_pair_name)
    {
        if(!(indices_vec.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] FlexibleLocalDihedralForceFieldGenerator : "
                   "parameter number of"
                   "indices_vec (" << indices_vec.size() << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string fld_expression =
            "c + ksin1*sin(  theta) + kcos1*cos(  theta)"
            "  + ksin2*sin(2*theta) + kcos2*cos(2*theta)"
            "  + ksin3*sin(3*theta) + kcos3*cos(3*theta)";
        auto torsion_ff = std::make_unique<OpenMM::CustomTorsionForce>(fld_expression);
        torsion_ff->addPerTorsionParameter("c");
        torsion_ff->addPerTorsionParameter("ksin1");
        torsion_ff->addPerTorsionParameter("kcos1");
        torsion_ff->addPerTorsionParameter("ksin2");
        torsion_ff->addPerTorsionParameter("kcos2");
        torsion_ff->addPerTorsionParameter("ksin3");
        torsion_ff->addPerTorsionParameter("kcos3");

        for(const auto& indices : indices_vec_)
        {
            torsion_ff->addTorsion(
                    indices[0], indices[1], indices[2], indices[3],
                    {fourier_table_[0],
                     fourier_table_[1], fourier_table_[2], fourier_table_[3],
                     fourier_table_[4], fourier_table_[5], fourier_table_[6]});
        }

        return torsion_ff;
    }

    void add_exclusion(std::vector<std::pair<std::size_t, std::size_t>>& exclusion_pairs) const noexcept
    {
        // TODO
        // duplication in exclusion list make error
        // all exclusion should be specified in HarnomicBond and GoContact

        // exclusion_pairs.push_back(std::make_pair(indices[0], indices[3]));
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::array<double, 7>     fourier_table_;
    std::string               aa_pair_name_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_DIHEDRAL_FORCE_FIELD_GENERATOR_HPP
