#ifndef OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>

#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"

#include "src/Topology.hpp"

HarmonicBondForceFieldGenerator
read_genesis_harmonic_bond_ff_generator(
        const std::vector<std::string>& bonds_data, Topology& topology, const bool use_periodic);

GaussianBondForceFieldGenerator
read_genesis_gaussian_bond_ff_generator(
        const std::vector<std::string>& angles_data, const bool use_periodic);

GoContactForceFieldGenerator
read_genesis_go_contact_ff_generator(
        const std::vector<std::string>& pairs_data, Topology& topology,
        const bool use_periodic);

FlexibleLocalAngleForceFieldGenerator
read_genesis_flexible_local_angle_ff_generator(
        const std::vector<std::string>& angles_data, const std::vector<std::string>& atoms_data,
        const std::string& aa_type, const bool use_periodic);

GaussianDihedralForceFieldGenerator
read_genesis_gaussian_dihedral_ff_generator(
        const std::vector<std::string>& dihedrals_data,
        const bool use_periodic);

FlexibleLocalDihedralForceFieldGenerator
read_genesis_flexible_local_dihedral_ff_generator(
        const std::vector<std::string>& dihedrals_data,
        const std::vector<std::string>& atoms_data,
        const std::pair<std::string, std::string>& aa_type_pair,
        const bool use_periodic);

ExcludedVolumeForceFieldGenerator
read_genesis_exv_ff_generator(const std::vector<std::string>& atomtypes_data,
        const std::vector<std::string>& atoms_data, Topology& topology,
        const bool use_periodic,
        const std::size_t ignore_particle_within_bond = 3);

#endif // OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
