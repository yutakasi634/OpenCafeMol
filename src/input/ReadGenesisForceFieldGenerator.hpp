#ifndef OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>

#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/CosineDihedralForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"

#include "src/Topology.hpp"

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_bonds_section(
        const std::vector<std::string>& bonds_data, Topology& topology,
        const bool use_periodic);

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_angles_section(
        const std::vector<std::string>& angles_data,
        const bool use_periodic, const std::vector<std::string>& res_name_vec);

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_dihedrals_section(
        const std::vector<std::string>& dihedrals_data,
        const bool use_periodic, const std::vector<std::string>& res_name_vec);

GoContactForceFieldGenerator
read_genesis_pairs_section(
        const std::vector<std::string>& pairs_data, Topology& topology,
        const bool use_periodic);

ExcludedVolumeForceFieldGenerator
read_genesis_exv_ff_generator(const std::vector<std::string>& atomtypes_data,
        const std::vector<std::string>& atoms_data, Topology& topology,
        const bool use_periodic,
        const std::size_t ignore_particle_within_bond = 3);

#endif // OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
