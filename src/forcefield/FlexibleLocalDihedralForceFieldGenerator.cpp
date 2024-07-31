#include "FlexibleLocalDihedralForceFieldGenerator.hpp"

#include "src/util/Constants.hpp"

std::unique_ptr<OpenMM::Force> FlexibleLocalDihedralForceFieldGenerator::generate() const
{
    double phi        = -Constant::pi;
    double min_energy = std::numeric_limits<double>::max();
    while(phi < Constant::pi)
    {
        const double sin1   = std::sin(phi);
        const double cos1   = std::cos(phi);
        const double sin_sq = sin1 * sin1;
        const double cos_sq = cos1 * cos1;
        const double sin2   = 2 * sin1 * cos1;
        const double cos2   = cos_sq - sin_sq;
        const double sin3   = sin1 * (3 - 4 * sin_sq);
        const double cos3   = cos1 * (4 * cos_sq - 3);

        min_energy = std::min(min_energy,
                        fourier_table_[1] * cos1 + fourier_table_[2] * sin1 +
                        fourier_table_[3] * cos2 + fourier_table_[4] * sin2 +
                        fourier_table_[5] * cos3 + fourier_table_[6] * sin3 +
                        fourier_table_[0]);
        phi += 1e-4;
    }

    const std::string fld_expression = fmt::format("{id}_k * ("
            "{id}_c - {id}_min_energy +"
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
    torsion_ff->addPerTorsionParameter(fmt::format("{}_min_energy", ffgen_id_));

    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const auto&  indices = indices_vec_[idx];
        const double k       = ks_[idx];
        torsion_ff->addTorsion(
                indices[0], indices[1], indices[2], indices[3],
                {k, fourier_table_[0],
                 fourier_table_[1], fourier_table_[2], fourier_table_[3],
                 fourier_table_[4], fourier_table_[5], fourier_table_[6],
                 min_energy});
    }

    return torsion_ff;
}

