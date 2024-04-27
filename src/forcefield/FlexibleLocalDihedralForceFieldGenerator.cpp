#include "FlexibleLocalDihedralForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> FlexibleLocalDihedralForceFieldGenerator::generate() const noexcept
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

