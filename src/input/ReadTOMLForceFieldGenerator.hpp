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
#include "src/forcefield/UniformLennardJonesAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/UniformWeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BasePairForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2CrossStackingForceFieldGenerator.hpp"
#include "src/forcefield/CombinatorialGoContactForceFieldGenerator.hpp"
#include "src/forcefield/PullingForceFieldGenerator.hpp"
#include "src/forcefield/PositionRestraintForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicCoMPullingForceFieldGenerator.hpp"

#include "src/Topology.hpp"
#include "src/toml11_fwd.hpp"

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

// TODO: enable to use offset
const ThreeSPN2BondForceFieldGenerator
read_toml_3spn2_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

const HarmonicAngleForceFieldGenerator
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

const CosineDihedralForceFieldGenerator
read_toml_cosine_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

const FlexibleLocalDihedralForceFieldGenerator
read_toml_flexible_local_dihedral_ff_generator(
        const toml::value& local_ff_data, const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table, Topology& topology,
        const bool use_periodic);

// TODO: enable to use offset
const ThreeSPN2BaseStackingForceFieldGenerator
read_toml_3spn2_base_stacking_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic);

// ----------------------------------------------------------------------------
// read global force field

const std::vector<std::pair<std::size_t, std::size_t>>
read_ignore_molecule_and_particles_within(const toml::value& ignore_table, const Topology& topology);

const std::vector<std::pair<std::string, std::string>>
read_ignore_group(const toml::value& ignore_table);

const ExcludedVolumeForceFieldGenerator
read_toml_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

const ThreeSPN2ExcludedVolumeForceFieldGenerator
read_toml_3spn2_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

// WCA input is like below
// [[forcefields.global]]
// potential   = "WCA"
// ignore.particles_within = {bond = 1, angle = 1}
// env.popc_epsilon = 0.416
// env.popc_sigma_T = 7.111
// env.popc_sigma_H = 4.62215
// parameters = [
//     {index =     0, sigma = "popc_sigma_H", epsilon = "popc_epsilon"},
//     # ...
// ]
const WeeksChandlerAndersenForceFieldGenerator
read_toml_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// UniformWeeksChandlerAndersen input is like below
//[[forcefields.global]]
//interaction = "Pair"
//potential   = "WCA"
//table.ASP.Head = {sigma =   4.1260, epsilon =   0.4996}
//# ...
//parameters = [
//    {index =     0, name = "Head"},
//    {index =     1, name =  "ASP"},
//    # ...
//]
const UniformWeeksChandlerAndersenForceFieldGenerator
read_toml_uniform_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>> ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// DebyeHuckel input is like below
// [[forcefields.global]]
// interaction = "Pair"
// potential   = "DebyeHuckel"
// ignore.particles_within.bond = 3
// parameters = [
//     {index =   2, charge = -0.6},
//     # ...
// ]
const DebyeHuckelForceFieldGenerator
read_toml_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

// iSoLFAttractive input is like below
// [[forcefields.global]]
// potential   = "iSoLFAttractive"
// ignore.particles_within = {bond = 1, angle = 1}
// env.popc_epsilon = 0.416
// env.popc_omega   = 9.867
// env.popc_sigma_T = 7.111
// parameters = [
//     {index =     2, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
//     # ...
// ]
const iSoLFAttractiveForceFieldGenerator
read_toml_isolf_attractive_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic);

// [[forcefields.global]]
// potential = "LennardJonesAttractive"
// cutoff    = 5.0
//table.A.A = {sigma =   4.1260, epsilon =  -0.7350}
//table.A.B = {sigma =   6.3361, epsilon =  -0.3710}
//table.B.B = {sigma =   4.1221, epsilon =   0.6500}
//parameters = [
//{index =     1, name = "A"},
//{index =     2, name = "A"},
//{index =     3, name = "B"},
//{index =     4, name = "B"},
// ...
// ]
const UniformLennardJonesAttractiveForceFieldGenerator
read_toml_uniform_lennard_jones_attractive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>> ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// TODO: enable to use offset
const ThreeSPN2BasePairLocalForceFieldGenerator
read_toml_3spn2_base_pair_local_ff_generator(
    const ThreeSPN2BasePairPotentialParameterBase& para,
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> base_pair,
    const bool use_periodic);

const ThreeSPN2BasePairForceFieldGenerator
read_toml_3spn2_base_pair_ff_generator(
        const ThreeSPN2BasePairPotentialParameterBase& para,
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> base_pair,
        const bool use_periodic);

// TODO: enable to use offset
const ThreeSPN2CrossStackingLocalForceFieldGenerator
read_toml_3spn2_cross_stacking_local_ff_generator(
    const ThreeSPN2CrossStackingPotentialParameterBase& para,
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
    const bool use_periodic
);

const ThreeSPN2CrossStackingForceFieldGenerator
read_toml_3spn2_cross_stacking_ff_generator(
        const ThreeSPN2CrossStackingPotentialParameterBase& para,
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
        const bool use_periodic);

// [[forcefields.global]]
// interaction = "Pair"
// potential   = "LennardJonesAttractive"
// cutoff      = 5.0
// parameters = [
// {index = 1, offset = 10, epsilon = 1.0, sigma = 1.0},
// ...
// ]
const LennardJonesAttractiveForceFieldGenerator
read_toml_lennard_jones_attractive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// [[forcefields.global]]
// interaction = "Pair"
// potential   = "LennardJonesRepulsive"
// cutoff      = 5.0
// parameters = [
// {index = 1, offset = 10, epsilon = 1.0, sigma = 1.0},
// ...
// ]
const LennardJonesRepulsiveForceFieldGenerator
read_toml_lennard_jones_repulsive_ff_generator(
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

// CombinatorialGocontact table is like below
// [[forcefields.global]]
// interaction = "CombinatorialGoContact"
// potential   = "CombinatorialGoContact"
// topology    = "contact"
// cutoff      = 2.0
// parameters = [
//     {indices_pair = [[ 0, 1, 2], [ 3, 4]], k = 0.3384, v0 = 9.3485},
//     # ...
// ]
// TODO: enable to add topology information
std::vector<CombinatorialGoContactForceFieldGenerator>
read_toml_combinatorial_go_contact_ff_generators(
        const toml::value& global_ff_data, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic);

// PDNS table is like below
// [[forcefields.global]]
// interaction = "PDNS"
// potential   = "PDNS"
// sigma = 1.0  # angstrom
// delta = 0.17453 # radian (= 10 degree)
// parameters  = [
// {index =    2, S3 = 1, kind = "DNA"},
// {index =    5, S3 = 4, kind = "DNA"},
// # ...
// {index = 1000, kind = "Protein", PN =  999, PC = 1001, k = 1.2, r0 = 5.0, theta0 = 100.0, phi0 = 130.0},
// {index = 1023, kind = "Protein", PN = 1022, PC = 1024, k = 1.2, r0 = 6.0, theta0 = 110.0, phi0 = 120.0},
// # ...
// ]
ProteinDNANonSpecificForceFieldGenerator
read_toml_protein_dna_non_specific_ff_generator(
    const  toml::value& global_ff_data, const bool use_periodic);

// -----------------------------------------------------------------------------
// read external force field

// PullingForce input is like below
// [[forcefields.external]]
// interaction = "PullingForce"
// parameters = [
//     {index = 0, force = [1.0, 0.0, 0.0]},
//     {index = 1, force = [0.0, 1.0, 0.0]},
//     # ...
// ]
PullingForceFieldGenerator
read_toml_pulling_ff_generator(
        const toml::value& external_ff_data, const Topology& topology,
        const bool use_periodic);


// PositionRestraint input is like below
// [[forcefields.external]]
// interaction = "PositionRestraint"
// potential   = "Harmonic"
// parameters = [
//     {index = 0, position = [0.0, 0.0, 0.0], k = 0.1, v0 = 10.0},
//     # ...
// ]
PositionRestraintForceFieldGenerator
read_toml_position_restraint_ff_generator(
        const toml::value& external_ff_data, const Topology& topology);

// HarmonicCoMPulling input is like below
// [[forcefields.external]]
// interaction = "CoMPulingForce"
// potential   = "Harmonic"
// parameters = [
//     {k = 10.0, v0 = 0.0, indices_pair = [[0, 1], [2, 3]]},
//     # ...
// ]
HarmonicCoMPullingForceFieldGenerator
read_toml_harmonic_com_pulling_ff_generator(
        const toml::value& external_ff_data, const bool use_periodic,
        const toml::value& env);

#endif // OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
