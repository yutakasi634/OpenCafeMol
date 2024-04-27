#ifndef OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include "Utility.hpp"

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
template<typename PotentialParameterType>
const ThreeSPN2BasePairLocalForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_base_pair_local_ff_generator(
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> base_pair,
    const bool use_periodic
)
{
    const std::string base0 = base_pair.first;
    const std::string base1 = base_pair.second;

    if (! ((base0 == "A" && base1 == "T") || (base0 == "T" && base1 == "A") ||
           (base0 == "G" && base1 == "C") || (base0 == "C" && base1 == "G")))
    {
        throw std::runtime_error(
            "[error] invalid base pair type " + base0 + "-" + base1 + " here. " +
            "One of the A-T, T-A, G-C, C-G is expected.");
    }

    // [[forcefields.local]]
    // interaction = "3SPN2BasePairLocal"
    // potential   = "3SPN2"
    // ignore.particles_within.nucleotide = 3
    // parameters  = [
    // {nucleotide_pair=[{       S = 0, B = 1, Base = "A"},{P = 187, S = 188, B = 189, Base = "T"}]},
    // {nucleotide_pair=[{P = 2, S = 3, B = 4, Base = "T"},{P = 184, S = 185, B = 186, Base = "A"}]},

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env    = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;

    for (const auto& param: params)
    {
        const auto nucleotide_pair = toml::find<toml::array>(param, "nucleotide_pair");

        if (nucleotide_pair.size() != 2)
        {
            throw std::runtime_error(
                "[error] invalid toml array size."
                "The array size of nucleotide_pair must be 2.");
        }

        const auto base_strand0 = toml::find<std::string>(nucleotide_pair[0], "Base");
        const auto base_strand1 = toml::find<std::string>(nucleotide_pair[1], "Base");

        if ((base_strand0 == base0) && (base_strand1 == base1))
        {
            indices_vec.push_back({
                Utility::find_parameter<size_t>(nucleotide_pair[0], env, "S"),
                Utility::find_parameter<size_t>(nucleotide_pair[0], env, "B"),
                Utility::find_parameter<size_t>(nucleotide_pair[1], env, "B"),
                Utility::find_parameter<size_t>(nucleotide_pair[1], env, "S")
            });
        }
        else if ((base_strand0 == base1) && (base_strand1 == base0))
        {
            indices_vec.push_back({
                Utility::find_parameter<size_t>(nucleotide_pair[1], env, "S"),
                Utility::find_parameter<size_t>(nucleotide_pair[1], env, "B"),
                Utility::find_parameter<size_t>(nucleotide_pair[0], env, "B"),
                Utility::find_parameter<size_t>(nucleotide_pair[0], env, "S")
            });
        }
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    3SPN2BasePairLocal   : " + PotentialParameterType::name + " BasePair "
          << base0 + "-" + base1
          << " (" << indices_vec.size() << " pairs found)"
          << std::endl;

    return ThreeSPN2BasePairLocalForceFieldGenerator<PotentialParameterType>(
        indices_vec, base_pair, use_periodic);
}

template<typename PotentialParameterType>
const ThreeSPN2BasePairForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_base_pair_ff_generator(
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> base_pair,
        const bool use_periodic)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "parameters", "env"});

    const std::string donor    = base_pair.first;
    const std::string acceptor = base_pair.second;

    if (! ((donor == "A" && acceptor == "T") || (donor == "T" && acceptor == "A") ||
           (donor == "G" && acceptor == "C") || (donor == "C" && acceptor == "G")))
    {
        throw std::runtime_error(
            "[error] invalid base pair type " + donor + "-" + acceptor + " here. " +
            "One of the A-T, T-A, G-C, C-G is expected.");
    }

    // [[forcefields.global]]
    // interaction = "3SPN2BasePair"
    // potential   = "3SPN2"
    // ignore.particles_within.nucleotide = 3
    // parameters = [
    //     {strand = 0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    const auto  pot    = toml::find<std::string>(global_ff_data, "potential");
    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env    = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    struct Nucleotide
    {
        static constexpr std::size_t nil() noexcept
        {
            return std::numeric_limits<std::size_t>::max();
        }

        Nucleotide() noexcept
            : strand(nil()), P(nil()), S(nil()), B(nil()), base("X")
        {}
        ~Nucleotide() noexcept = default;
        Nucleotide(Nucleotide const&) = default;
        Nucleotide(Nucleotide &&)     = default;
        Nucleotide& operator=(Nucleotide const&) = default;
        Nucleotide& operator=(Nucleotide &&)     = default;

        std::size_t strand;
        std::size_t P, S, B;
        std::string base;
    };

    // parameters = [
    //     {strand = 0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    std::vector<Nucleotide> nucleotide_idxs;
    for(const auto& param : params)
    {
        Nucleotide nucleotide;

        if(param.as_table().count("P") != 0)
        {
            nucleotide.P  = Utility::find_parameter<size_t>(param, env, "P");
        }
        nucleotide.S      = Utility::find_parameter<size_t>(param, env, "S");
        nucleotide.B      = Utility::find_parameter<size_t>(param, env, "B");
        nucleotide.strand = Utility::find_parameter<size_t>(param, env, "strand");

        const auto base = toml::find<std::string>(param, "Base");
        if (base=="A" || base == "T" || base == "G" || base == "C")
        {
            nucleotide.base = base;
        }
        else
        {
            throw std::runtime_error(
                "[error] invalid base type " + base + " here. "
                "One of the \"A\", \"T\", \"C\", \"G\" is expected."
            );
        }

        nucleotide_idxs.push_back(nucleotide);
    }

    // Create base-pairing list
    std::vector<std::array<std::size_t, 2>> indices_donor;     // 0: base, 1: sugar
    std::vector<std::array<std::size_t, 2>> indices_acceptor;  // 0: base, 1: sugar

    for(const auto& nuc : nucleotide_idxs)
    {
        if (nuc.base == donor)
        {
            indices_donor.push_back({nuc.B, nuc.S});
        }
        else if (nuc.base == acceptor)
        {
            indices_acceptor.push_back({nuc.B, nuc.S});
        }
    }

    // ignore list generation
    index_pairs_type ignore_list;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list = read_ignore_molecule_and_particles_within(ignore, topology);
    }

    std::cerr << "    Global        : " + PotentialParameterType::name + " BasePair "
          << donor << "-" << acceptor << " ("
          << indices_donor.size()    << donor << " and "
          << indices_acceptor.size() << acceptor
          << " found)" << std::endl;

    return ThreeSPN2BasePairForceFieldGenerator<PotentialParameterType>(
        indices_donor, indices_acceptor, base_pair, ignore_list,
        use_periodic);
}



// TODO: enable to use offset
template<typename PotentialParameterType>
const ThreeSPN2CrossStackingLocalForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_cross_stacking_local_ff_generator(
    const toml::value& local_ff_data, Topology& topology,
    const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
    const bool use_periodic
)
{
    // ================================================================
    // cross stacking
    //
    //  Sense                   Anti-sense
    //    5'                        3'
    //    ^    /               \    |
    //    | P o                o P  |  P : a phosphate site
    //    |    \    B0    Bp   /    |  S : a sugar site
    //    |  S  o -- o===o -- o S   |  B0: the reference nucleobase
    //    |    /      \ /     \     |  Bp: the base-pairing nucleobase
    //    | P o        X       o P  |
    //    |    \      / \     /     |
    //    |     o -- o===o -- o     |  Bn: the neighboring nucleobase
    //    |    /    Bn   Bc   \     |  Bc: the cross-stacking nucleobase
    //    | P o                o P  |
    //    |    \              /     v
    //    3'                        5'

    // [[forcefields.local]]
    // interaction = "3SPN2CrossStackingLocal"
    // potential   = "3SPN2"
    // ignore.particles_within.nucleotide = 3
    // parameters  = [
    // {nucleotide_group=[{P =   2, S =   3, B =   4, Base = "T"}, {P = 184, S = 185, B = 186, Base = "A"},
    //                    {P =   5, S =   6, B =   7, Base = "A"}, {P = 187, S = 188, B = 189, Base = "T"}]},
    // here, the order of nucleotide_group must be
    //  {nucleotide_group = [{the reference   nucleotide}, {the base pairing   nucleotide},
    //                       {the neighboring nucleotide}, {the cross-stacking nucleotide}]}

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env    = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 6>> indices_vec;
    std::vector<std::string>                base_kind_vec;

    for (const auto& param: params)
    {
        const auto nucleotide_group = toml::find<toml::array>(param, "nucleotide_group");

        if (nucleotide_group.size() != 4)
        {
            throw std::runtime_error(
                "[error] invalid toml array size."
                "The array size of nucleotide_group must be 4.");
        }

        const auto base_0 = toml::find<std::string>(nucleotide_group[0], "Base"); // reference
        const auto base_p = toml::find<std::string>(nucleotide_group[1], "Base"); // base-pairing
        const auto base_n = toml::find<std::string>(nucleotide_group[2], "Base"); // neighboring
        const auto base_c = toml::find<std::string>(nucleotide_group[3], "Base"); // cross-stacking

        if(strand_kind == "sense") // Base^{5'}↑
        {
            if (base_0 == bp_kind.first && base_p == bp_kind.second)
            {
                indices_vec.push_back({                                // d: donor, a:acceptor
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "S"), // p1 (=d2)
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "B"), // p2 (=d1)
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "B"), // p3 (=a1)
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "S"), // p4 (=a2)
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "B"), // p5 (=d3)
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "B"), // p6 (=a3)
                });
                base_kind_vec.push_back(base_0 + base_c + "5");
            }
            if (base_c == bp_kind.first && base_n == bp_kind.second)
            {
                indices_vec.push_back({                                // d: donor, a:acceptor
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "S"), // p1 (=d2)
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "B"), // p2 (=d1)
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "B"), // p3 (=a1)
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "S"), // p4 (=a2)
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "B"), // p5 (=d3)
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "B"), // p6 (=a3)
                });
                base_kind_vec.push_back(base_c + base_0 + "5");
            }
        }
        else if (strand_kind == "antisense") // Base^{3'}↓
        {
            if (base_0 == bp_kind.first && base_p == bp_kind.second)
            {
                indices_vec.push_back({                                // d: donor, a:acceptor
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "S"), // p1 (=d2)
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "B"), // p2 (=d1)
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "B"), // p3 (=a1)
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "S"), // p4 (=a2)
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "B"), // p5 (=d3)
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "B"), // p6 (=a3)
                });
                base_kind_vec.push_back(base_p + base_n + "3");
            }
            if (base_c == bp_kind.first && base_n == bp_kind.second)
            {
                indices_vec.push_back({                                // d: donor, a:acceptor
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "S"), // p1 (=d2)
                    Utility::find_parameter<size_t>(nucleotide_group[2], env, "B"), // p2 (=d1)
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "B"), // p3 (=a1)
                    Utility::find_parameter<size_t>(nucleotide_group[3], env, "S"), // p4 (=a2)
                    Utility::find_parameter<size_t>(nucleotide_group[0], env, "B"), // p5 (=d3)
                    Utility::find_parameter<size_t>(nucleotide_group[1], env, "B"), // p6 (=a3)
                });
                base_kind_vec.push_back(base_n + base_p + "3");
            }
        }
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    3SPN2CrossStackingLocal   : " + PotentialParameterType::name + " CrossStacking "
              << bp_kind.first + "-" + bp_kind.second <<" in "
              << strand_kind << " strand "
              << " (" << indices_vec.size() << " pairs found)"
              << std::endl;

    return ThreeSPN2CrossStackingLocalForceFieldGenerator<PotentialParameterType>(
        indices_vec, base_kind_vec, bp_kind, use_periodic);
}

template<typename PotentialParameterType>
const ThreeSPN2CrossStackingForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_cross_stacking_ff_generator(
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
        const bool use_periodic)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "parameters", "env"});

    if (! ((bp_kind.first == "A" && bp_kind.second == "T" )||
           (bp_kind.first == "T" && bp_kind.second == "A" )||
           (bp_kind.first == "G" && bp_kind.second == "C" )||
           (bp_kind.first == "C" && bp_kind.second == "G" )))
    {
        throw std::runtime_error(
            "[error] invalid base type " + bp_kind.first + '-' + bp_kind.second + " here. " +
            "One of the A-T, T-A, G-C, C-G is expected."
        );
    }

    // [[forcefields.global]]
    // interaction = "3SPN2CrossStacking"
    // potential   = "3SPN2"
    // ignore.particles_within.nucleotide = 3
    // parameters = [
    //     {strand = 0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    const auto  pot    = toml::find<std::string>(global_ff_data, "potential");
    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env    = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    struct Nucleotide
    {
        static constexpr std::size_t nil() noexcept
        {
            return std::numeric_limits<std::size_t>::max();
        }

        Nucleotide() noexcept
            : strand(nil()), P(nil()), S(nil()), B(nil()), base("X")
        {}
        ~Nucleotide() noexcept = default;
        Nucleotide(Nucleotide const&) = default;
        Nucleotide(Nucleotide &&)     = default;
        Nucleotide& operator=(Nucleotide const&) = default;
        Nucleotide& operator=(Nucleotide &&)     = default;

        std::size_t strand;
        std::size_t P, S, B;
        std::string base;
    };

    // parameters = [
    //     {strand = 0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    std::vector<Nucleotide> nucleotide_idxs;
    for(const auto& param: params)
    {
        Nucleotide nucleotide;

        if(param.as_table().count("P") != 0)
        {
            nucleotide.P  = Utility::find_parameter<size_t>(param, env, "P");
        }
        nucleotide.S      = Utility::find_parameter<size_t>(param, env, "S");
        nucleotide.B      = Utility::find_parameter<size_t>(param, env, "B");
        nucleotide.strand = Utility::find_parameter<size_t>(param, env, "strand");

        const auto base = toml::find<std::string>(param, "Base");
        if (base=="A" || base == "T" || base == "G" || base == "C")
        {
            nucleotide.base = base;
        }
        else
        {
            throw std::runtime_error(
                "[error] invalid base type " + base + " here. "
                "One of the \"A\", \"T\", \"C\", \"G\" is expected."
            );
        }

        nucleotide_idxs.push_back(nucleotide);
    }

    // ================================================================
    // cross stacking
    //
    //  Sense                   Anti-sense
    //    5'                        3'
    //    ^    /               \    |
    //    | P o                o P  |  P : a phosphate site
    //    |    \    B0    Bp   /    |  S : a sugar site
    //    |  S  o -- o===o -- o S   |  B0: the reference nucleobase
    //    |    /      \ /     \     |  Bp: the base-pairing nucleobase
    //    | P o        X       o P  |
    //    |    \      / \     /     |
    //    |     o -- o===o -- o     |  Bn: the neighboring nucleobase
    //    |    /    Bn   Bc   \     |  Bc: the cross-stacking nucleobase
    //    | P o                o P  |
    //    |    \              /     v
    //    3'                        5'

    // Create cross-stacking base-pair list
    std::vector<std::array<std::size_t, 3>> indices_donor;    // {B0, S, Bn}
    std::vector<std::array<std::size_t, 3>> indices_acceptor; // {Bp, S, Bc}
    std::vector<std::string>                base_kind_acceptor;   // {Bp, Bc, "(5|3)"}

    if(strand_kind == "sense") // Base^{5'}↑
    {
        for (std::size_t i = 0; i < nucleotide_idxs.size(); ++i)
        {
            const auto &nuc = nucleotide_idxs.at(i);

            // Donor: reference base (B0) and neighboring base (Bn)
            if (nuc.base == bp_kind.first)
            {
                if (i + 1 < nucleotide_idxs.size() &&
                    nucleotide_idxs.at(i + 1).strand == nuc.strand)
                {
                    const auto &nuc3 = nucleotide_idxs.at(i + 1);
                    indices_donor.push_back({nuc.B, nuc.S, nuc3.B});
                }
            }
            // Acceptor: base pair (Bp) and crossing base (Bc)
            if (nuc.base == bp_kind.second)
            {
                if (0 < i && nucleotide_idxs.at(i - 1).strand == nuc.strand)
                {
                    const auto &nuc5 = nucleotide_idxs.at(i - 1);
                    indices_acceptor.push_back({nuc.B, nuc.S, nuc5.B});
                    base_kind_acceptor.push_back(bp_kind.first + nuc5.base + "5");
                }
            }
        }
    }
    else if (strand_kind == "antisense") // Base^{3'}↓
    {
        for (std::size_t i = 0; i < nucleotide_idxs.size(); ++i)
        {
            const auto &nuc = nucleotide_idxs.at(i);

            // Acceptor: base pair (Bp) and crossing base (Bc)
            if (nuc.base == bp_kind.first)
            {
                if (i + 1 < nucleotide_idxs.size() &&
                    nucleotide_idxs.at(i + 1).strand == nuc.strand)
                {
                    const auto &nuc3 = nucleotide_idxs.at(i + 1);
                    indices_acceptor.push_back({nuc.B, nuc.S, nuc3.B});
                    base_kind_acceptor.push_back(bp_kind.second  + nuc3.base + "3");
                }
            }
            // Donor: reference base (B0) and neighboring base (Bn)
            if (nuc.base == bp_kind.second)
            {
                if (0 < i && nucleotide_idxs.at(i - 1).strand == nuc.strand)
                {
                    const auto &nuc5 = nucleotide_idxs.at(i - 1);
                    indices_donor.push_back({nuc.B, nuc.S, nuc5.B});
                }
            }
        }
    }
    else
    {
        throw std::runtime_error(
            "[error] invalid strand type " + strand_kind + " here. " +
            "\"sense\" or \"antisense\" is expected."
        );
    }

    // ignore list generation
    index_pairs_type ignore_list;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list = read_ignore_molecule_and_particles_within(ignore, topology);
    }

    std::cerr << "    Global        : " + PotentialParameterType::name + " CrossStacking "
              << bp_kind.first << "-" << bp_kind.second << " in "
              << strand_kind << " strand ("
              << indices_donor.size()    << " x "
              << indices_acceptor.size() << " pairs found)" << std::endl;

    return ThreeSPN2CrossStackingForceFieldGenerator<PotentialParameterType>(
        indices_donor, indices_acceptor, base_kind_acceptor, bp_kind, ignore_list,
        use_periodic);
}

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
