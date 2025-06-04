#ifndef OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BondForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicAngleForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "src/forcefield/CosineDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BaseStackingForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BasePairLocalForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2CrossStackingLocalForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/WeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/DebyeHuckelForceFieldGenerator.hpp"
#include "src/forcefield/iSoLFAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/LennardJonesAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/LennardJonesRepulsiveForceFieldGenerator.hpp"
#include "src/forcefield/TrigonometricForceFieldGenerator.hpp"
#include "src/forcefield/UniformLennardJonesAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/UniformWeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BasePairForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2CrossStackingForceFieldGenerator.hpp"
#include "src/forcefield/CombinatorialGoContactForceFieldGenerator.hpp"
#include "src/forcefield/ProteinDNANonSpecificForceFieldGenerator.hpp"
#include "src/forcefield/PullingForceFieldGenerator.hpp"
#include "src/forcefield/PositionRestraintForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicCoMPullingForceFieldGenerator.hpp"
#include "src/forcefield/EXVRectangularBoxForceFieldGenerator.hpp"
#include "src/forcefield/CappedGoContactForceFieldGenerator.hpp"

#include "src/Topology.hpp"
#include "src/input/Utility.hpp"
#include "src/util/Utility.hpp"

#include <OpenMM.h>

// -----------------------------------------------------------------------------
// read local force field

HarmonicBondForceFieldGenerator
read_toml_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic);

GaussianBondForceFieldGenerator
read_toml_gaussian_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

GoContactForceFieldGenerator
read_toml_go_contact_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

CappedGoContactForceFieldGenerator
read_toml_capped_go_contact_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

// TODO: enable to use offset
ThreeSPN2BondForceFieldGenerator
read_toml_3spn2_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

HarmonicAngleForceFieldGenerator
read_toml_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic);

FlexibleLocalAngleForceFieldGenerator
read_toml_flexible_local_angle_ff_generator(
        const toml::value& local_ff_data, const std::string& aa_type,
        const std::array<double, 10> spline_table, Topology& topology,
        const bool use_periodic);

GaussianDihedralForceFieldGenerator
read_toml_gaussian_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

CosineDihedralForceFieldGenerator
read_toml_cosine_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

FlexibleLocalDihedralForceFieldGenerator
read_toml_flexible_local_dihedral_ff_generator(
        const toml::value& local_ff_data, const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table, Topology& topology,
        const bool use_periodic);

// TODO: enable to use offset
ThreeSPN2BaseStackingForceFieldGenerator
read_toml_3spn2_base_stacking_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

// ----------------------------------------------------------------------------
// read global force field

std::vector<std::pair<std::size_t, std::size_t>>
read_ignore_molecule_and_particles_within(const toml::value& ignore_table, const Topology& topology);

std::vector<std::pair<std::string, std::string>>
read_ignore_group(const toml::value& ignore_table);

ExcludedVolumeForceFieldGenerator
read_toml_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

ThreeSPN2ExcludedVolumeForceFieldGenerator
read_toml_3spn2_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

WeeksChandlerAndersenForceFieldGenerator
read_toml_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

UniformWeeksChandlerAndersenForceFieldGenerator
read_toml_uniform_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>> ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

DebyeHuckelForceFieldGenerator
read_toml_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

iSoLFAttractiveForceFieldGenerator
read_toml_isolf_attractive_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

UniformLennardJonesAttractiveForceFieldGenerator
read_toml_uniform_lennard_jones_attractive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>> ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// TODO: enable to use offset
ThreeSPN2BasePairLocalForceFieldGenerator
read_toml_3spn2_base_pair_local_ff_generator(
    const ThreeSPN2BasePairPotentialParameterBase& para,
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> base_pair,
    const bool use_periodic);

ThreeSPN2BasePairForceFieldGenerator
read_toml_3spn2_base_pair_ff_generator(
        const ThreeSPN2BasePairPotentialParameterBase& para,
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> base_pair,
        const bool use_periodic);

// TODO: enable to use offset
ThreeSPN2CrossStackingLocalForceFieldGenerator
read_toml_3spn2_cross_stacking_local_ff_generator(
    const ThreeSPN2CrossStackingPotentialParameterBase& para,
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
    const bool use_periodic
);

ThreeSPN2CrossStackingForceFieldGenerator
read_toml_3spn2_cross_stacking_ff_generator(
        const ThreeSPN2CrossStackingPotentialParameterBase& para,
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
        const bool use_periodic);

LennardJonesAttractiveForceFieldGenerator
read_toml_lennard_jones_attractive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

LennardJonesRepulsiveForceFieldGenerator
read_toml_lennard_jones_repulsive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

TrigonometricForceFieldGenerator
read_toml_trigonometric_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

CombinatorialGoContactForceFieldGenerator
read_toml_combinatorial_go_contact_ff_generator(
        const double cutoff_ratio, const toml::value& contacts_info,
        const toml::value& env, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>>& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>>& ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec);

std::vector<CombinatorialGoContactForceFieldGenerator>
read_toml_combinatorial_go_contact_ff_generators(
        const toml::value& global_ff_data, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

ProteinDNANonSpecificForceFieldGenerator
read_toml_protein_dna_non_specific_ff_generator(
    const  toml::value& global_ff_data, const bool use_periodic);

// -----------------------------------------------------------------------------
// read external force field
PullingForceFieldGenerator
read_toml_pulling_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic);

PositionRestraintForceFieldGenerator
read_toml_position_restraint_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic);

HarmonicCoMPullingForceFieldGenerator
read_toml_harmonic_com_pulling_ff_generator(
        const toml::value& external_ff_data, const bool use_periodic,
        const toml::value& env);

EXVRectangularBoxForceFieldGenerator
read_toml_exv_rectangular_box_ff_generator(
        const toml::value& external_ff_data, const bool use_periodic,
        const toml::value& env);

#endif // OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
