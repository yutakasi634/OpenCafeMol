#ifndef OPEN_AICG2_PLUS_READ_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/DebyeHuckelForceFieldGenerator.hpp"

const HarmonicBondForceFieldGenerator
read_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    for(const auto& param : params)
    {
        const auto&  indices =
            toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
        const double v0 =
            toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
        const double k =
            toml::find<double>(param, "k") * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm; // KJ/(mol nm^2)

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }
    topology.add_edges(indices_vec, "bond");

    std::cerr << "    BondLength    : Harmonic" << std::endl;
    return HarmonicBondForceFieldGenerator(indices_vec, v0s, ks);
}

const GaussianBondForceFieldGenerator
read_gaussian_bond_ff_generator(const toml::value& local_ff_data)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;
    std::vector<double>                              sigmas;

    for(const auto& param : params)
    {
        const auto& indices =
            toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
        const double k  =
            toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double v0 =
            toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
        const double sigma =
            toml::get<double>(
                    find_either(param, "sigma", "σ")) * OpenMM::NmPerAngstrom; // nm

        indices_vec.push_back(indices);
        ks         .push_back(k);
        v0s        .push_back(v0);
        sigmas     .push_back(sigma);
    }

    std::cerr << "    BondLength    : Gaussian" << std::endl;
    return GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas);
}

const GoContactForceFieldGenerator
read_go_contact_ff_generator(const toml::value& local_ff_data, Topology& topology)
{
    // TODO: enable to optimization based on cutoff
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(const auto& param : params)
    {
        const auto& indices =
            toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
        const double k =
            toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double r0 =
            toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm

        ks         .push_back(k);
        indices_vec.push_back(indices);
        r0s        .push_back(r0);
    }
    topology.add_edges(indices_vec, "contact");

    std::cerr << "    BondLength    : GoContact" << std::endl;
    return GoContactForceFieldGenerator(indices_vec, ks, r0s);
}

const FlexibleLocalAngleForceFieldGenerator
read_flexible_local_angle_ff_generator(const toml::value& local_ff_data,
        const std::string& aa_type, const std::array<double, 10> spline_table)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::vector<std::size_t>> indices_vec;
    std::vector<double>                   ks;

    for(const auto& param : params)
    {
        const std::string y = toml::find<std::string>(param, "y");
        if(y.substr(3, 3) == aa_type) // y is like "y1_PHE"
        {
            const auto& indices =
                toml::find<std::vector<std::size_t>>(param, "indices");
            const double k =
                toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol

            indices_vec.push_back(indices);
            ks         .push_back(k);
        }
    }

    std::cerr << "    BondAngle     : FlexibleLocalAngle ("
              << aa_type << ")" << std::endl;

    return FlexibleLocalAngleForceFieldGenerator(
               indices_vec, ks, spline_table, aa_type);
}

const GaussianDihedralForceFieldGenerator
read_gaussian_dihedral_ff_generator(const toml::value& local_ff_data)
{

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     sigmas;

    for(const auto& param : params)
    {
        const auto& indices =
            toml::find<std::array<std::size_t, 4>>(param, "indices");
        const double k  =
            toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0 = toml::find<double>(param, "v0"); // radiuns
        const double sigma =
            toml::get<double>(find_either(param, "sigma", "σ")); // radiuns

        indices_vec.push_back(indices);
        ks         .push_back(k);
        theta0s    .push_back(theta0);
        sigmas     .push_back(sigma);
    }

    std::cerr << "    DihedralAngle : Gaussian" << std::endl;
    return GaussianDihedralForceFieldGenerator(indices_vec, ks, theta0s, sigmas);
}

const FlexibleLocalDihedralForceFieldGenerator
read_flexible_local_dihedral_ff_generator(const toml::value& local_ff_data,
        const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::string                             aa_pair_name;

    if(aa_pair_type.second == "GLY") // R1-GLY case
    {
        for(const auto& param : params)
        {
            const std::string coef    = toml::find<std::string>(param, "coef");
            const auto&       indices =
                toml::find<std::array<std::size_t, 4>>(param, "indices");
            const double      k       = toml::find<double>(param, "k");

            if(coef.substr(4, 3) == "GLY")
            {
                indices_vec.push_back(indices);
                ks         .push_back(k);
            }
        }
        aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
    }
    else if(aa_pair_type.second == "PRO")
    {
        if(aa_pair_type.first == "GLY") // GLY-PRO case
        {
            for(const auto& param : params)
            {
                const std::string coef    =
                    toml::find<std::string>(param, "coef");
                const auto&       indices =
                    toml::find<std::array<std::size_t, 4>>(param, "indices");
                const double      k       = toml::find<double>(param, "k");

                if(coef == "GLY-PRO")
                {
                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                }
            }
            aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
        }
        else
        {
            for(const auto& param : params)
            {
                const std::string coef    =
                    toml::find<std::string>(param, "coef");
                const auto&       indices =
                    toml::find<std::array<std::size_t, 4>>(param, "indices");
                const double      k       = toml::find<double>(param, "k");

                if(coef.substr(4, 3) == "PRO") // R2-PRO case
                {
                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                }
            }
            aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
        }
    }
    else // R1-R3 case
    {
        for(const auto& param : params)
        {
            const std::string coef    =
                toml::find<std::string>(param, "coef");
            const auto&       indices =
                toml::find<std::array<std::size_t, 4>>(param, "indices");
            const double      k       = toml::find<double>(param, "k");

            const std::string second_aa = coef.substr(4, 3);
            if(second_aa != "GLY" && second_aa != "PRO")
            {
                if(coef.substr(0, 3) == aa_pair_type.first)
                {
                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                }
            }
        }
        aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
    }

    std::cerr << "    DihedralAngle : FlexibleLocalDihedral ("
              << aa_pair_name << ")" << std::endl;
    return FlexibleLocalDihedralForceFieldGenerator(
               indices_vec, ks, fourier_table, aa_pair_name);
}

const ExcludedVolumeForceFieldGenerator
read_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology)
{
    const double eps =
        toml::find<double>(global_ff_data, "epsilon") * OpenMM::KJPerKcal; // KJPermol
    const double cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    const auto&  params = toml::find<toml::array>(global_ff_data, "parameters");
    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index  = toml::find<std::size_t>(param, "index");
        const double      radius =
            toml::find<double>(param, "radius") * OpenMM::NmPerAngstrom; // nm
        radius_vec.at(index) = radius;
    }

    ExcludedVolumeForceFieldGenerator::index_pairs_type ignore_list =
        topology.ignore_list_within_edge(3, "bond");
    ExcludedVolumeForceFieldGenerator::index_pairs_type contact_ignore_list =
        topology.ignore_list_within_edge(1, "contact");
    ignore_list.insert(ignore_list.end(),
          contact_ignore_list.begin(), contact_ignore_list.end());

    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    std::cerr << "    Global        : ExcludedVolume" << std::endl;
    return ExcludedVolumeForceFieldGenerator(eps, cutoff, radius_vec, ignore_list);
}

const DebyeHuckelForceFieldGenerator
read_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const DebyeHuckelForceFieldGenerator::index_pairs_type& contacted_pairs)
{
    using index_pairs_type = DebyeHuckelForceFieldGenerator::index_pairs_type;

    const double cutoff = toml::find_or(global_ff_data, "cutoff", 5.5);

    bool        ignore_molecule_flag = false;
    index_pairs_type ignore_list;

    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        if(ignore.contains("molecule"))
        {
            const std::string name = toml::find<std::string>(ignore, "molecule");
            if(name == "Intra" || name == "Self")
            {
                ignore_molecule_flag = true;
                index_pairs_type mol_ignore_list =
                    topology.ignore_list_within_molecule();
                ignore_list.insert(ignore_list.end(),
                     mol_ignore_list.begin(), mol_ignore_list.end());
            }
            else if(name == "Others" || name == "Inter")
            {
                throw std::runtime_error(
                    "[error] ignore molecule for the DebyeHuckelPotential do not support "
                    "\"Others\" or \"Inter\".");
            }
        }

        if(ignore.contains("particle_within"))
        {
            const auto particle_within = toml::find(ignore, "particle_within");
            if(particle_within.contains("bond"))
            {
                if(ignore_molecule_flag)
                {
                    std::cerr <<
                        "[warning] ignore molecule \"Intra\" or \"Self\" was defined,"
                        "so this ignore particle within bond will be ignored." << std::endl;
                }
                else
                {
                    std::size_t ignore_bond_within =
                        toml::find<std::size_t>(particle_within, "bond");
                    index_pairs_type bond_ignore_list =
                        topology.ignore_list_within_edge(ignore_bond_within, "bond");
                    ignore_list.insert(ignore_list.end(),
                        bond_ignore_list.begin(), bond_ignore_list.end());
                }
            }
            if(particle_within.contains("contact"))
            {
                std::size_t ignore_contact_within =
                    toml::find<std::size_t>(particle_within, "contact");
                index_pairs_type contact_ignore_list =
                    topology.ignore_list_within_edge(ignore_contact_within, "contact");
                ignore_list.insert(ignore_list.end(),
                    contact_ignore_list.begin(), contact_ignore_list.end());
            }
        }
    }
    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    std::vector<std::optional<double>> charge_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index  = toml::find<std::size_t>(param, "index");
        const double      charge = toml::find<double>(param, "charge");
        charge_vec.at(index) = charge;
    }

    std::cerr << "    Global        : DebyeHuckel" << std::endl;
    return DebyeHuckelForceFieldGenerator(
            ionic_strength, temperature, cutoff, charge_vec, ignore_list);
}

#endif // OPEN_AICG2_PLUS_READ_FORCE_FIELD_GENERATOR_HPP
