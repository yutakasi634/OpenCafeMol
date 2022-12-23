#ifndef OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicAngleForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/WeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/DebyeHuckelForceFieldGenerator.hpp"
#include "src/forcefield/iSoLFAttractiveForceFieldGenerator.hpp"

// -----------------------------------------------------------------------------
// read local force field

const HarmonicBondForceFieldGenerator
read_toml_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    for(const auto& param : params)
    {
        const auto&  indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const double v0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // nm
        const double k =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm; // KJ/(mol nm^2)

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : Harmonic (" << indices_vec.size() << " found)" << std::endl;
    return HarmonicBondForceFieldGenerator(indices_vec, v0s, ks);
}

const GaussianBondForceFieldGenerator
read_toml_gaussian_bond_ff_generator(const toml::value& local_ff_data, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;
    std::vector<double>                              sigmas;

    for(const auto& param : params)
    {
        const auto& indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const double k  =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double v0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // nm
        const double sigma =
            Utility::find_parameter_either<double>(param, env, "sigma", "σ")
            * OpenMM::NmPerAngstrom; // nm

        indices_vec.push_back(indices);
        ks         .push_back(k);
        v0s        .push_back(v0);
        sigmas     .push_back(sigma);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : Gaussian (" << indices_vec.size() << " found)" << std::endl;
    return GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas);
}

const GoContactForceFieldGenerator
read_toml_go_contact_ff_generator(const toml::value& local_ff_data, Topology& topology)
{
    // TODO: enable to optimization based on cutoff
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(const auto& param : params)
    {
        const auto& indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const double k =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double r0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // nm

        ks         .push_back(k);
        indices_vec.push_back(indices);
        r0s        .push_back(r0);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : GoContact (" << indices_vec.size() << " found)" << std::endl;
    return GoContactForceFieldGenerator(indices_vec, ks, r0s);
}

const HarmonicAngleForceFieldGenerator
read_toml_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     v0s;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        const auto&  indices =
            Utility::find_parameter<std::array<std::size_t, 3>>(param, env, "indices");
        double v0 =
            Utility::find_parameter<double>(param, env, "v0"); // radian
        if(v0 < 0 || Constant::pi*2 < v0)
        {
            throw std::runtime_error("[error] read_toml_harmonic_angle_ff_generator: "
                "v0 must be between 0 and 2pi");
        }
        else if(Constant::pi < v0)
        {
            v0 = 2.0 * Constant::pi - v0;
        }
        const double k =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondAngle     : Harmonic (" << indices_vec.size() << " found)" << std::endl;
    return HarmonicAngleForceFieldGenerator(indices_vec, v0s, ks);

}

const FlexibleLocalAngleForceFieldGenerator
read_toml_flexible_local_angle_ff_generator(const toml::value& local_ff_data,
        const std::string& aa_type, const std::array<double, 10> spline_table, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                   ks;

    for(const auto& param : params)
    {
        const std::string y = toml::find<std::string>(param, "y");
        if(y.substr(3, 3) == aa_type) // y is like "y1_PHE"
        {
            const auto& indices =
                Utility::find_parameter<std::array<std::size_t, 3>>(param, env, "indices");
            const double k =
                Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol

            indices_vec.push_back(indices);
            ks         .push_back(k);
        }
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondAngle     : FlexibleLocalAngle - "
              << aa_type << " (" << indices_vec.size() << " found)" << std::endl;

    return FlexibleLocalAngleForceFieldGenerator(
               indices_vec, ks, spline_table, aa_type);
}

const GaussianDihedralForceFieldGenerator
read_toml_gaussian_dihedral_ff_generator(const toml::value& local_ff_data, Topology& topology)
{

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     sigmas;

    for(const auto& param : params)
    {
        const auto& indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
        const double k  =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0 =
            Utility::find_parameter<double>(param, env, "v0"); // radiuns
        const double sigma =
            Utility::find_parameter_either<double>(param, env, "sigma", "σ"); // radiuns

        indices_vec.push_back(indices);
        ks         .push_back(k);
        theta0s    .push_back(theta0);
        sigmas     .push_back(sigma);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    DihedralAngle : Gaussian (" << indices_vec.size() << " found)" << std::endl;
    return GaussianDihedralForceFieldGenerator(indices_vec, ks, theta0s, sigmas);
}

const FlexibleLocalDihedralForceFieldGenerator
read_toml_flexible_local_dihedral_ff_generator(const toml::value& local_ff_data,
        const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table, Topology& topology)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

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
            const double      k       = Utility::find_parameter<double>(param, env, "k");

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
                const std::string coef    = toml::find<std::string>(param, "coef");
                const auto&       indices =
                    Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
                const double      k       =
                    Utility::find_parameter<double>(param, env, "k");

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
                const std::string coef    = toml::find<std::string>(param, "coef");
                const auto&       indices =
                    Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
                const double      k       =
                    Utility::find_parameter<double>(param, env, "k");

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
            const std::string coef    = toml::find<std::string>(param, "coef");
            const auto&       indices =
                Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
            const double      k       =
                Utility::find_parameter<double>(param, env, "k");

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

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    DihedralAngle : FlexibleLocalDihedral - "
              << aa_pair_name << " (" << indices_vec.size() << " found)" << std::endl;
    return FlexibleLocalDihedralForceFieldGenerator(
               indices_vec, ks, fourier_table, aa_pair_name);
}

// ----------------------------------------------------------------------------
// read global force field

const std::vector<std::pair<std::size_t, std::size_t>>
read_ignore_molecule_and_particles_within(const toml::value& ignore_table, const Topology& topology)
{
    std::vector<std::pair<std::size_t, std::size_t>> ignore_list;
    bool                                             ignore_molecule_flag = false;

    if(ignore_table.contains("molecule"))
    {
        const std::string name = toml::find<std::string>(ignore_table, "molecule");
        if(name == "Intra" || name == "Self")
        {
            ignore_molecule_flag = true;
            std::vector<std::pair<std::size_t, std::size_t>> mol_ignore_list =
                topology.ignore_list_within_molecule();
            ignore_list.insert(ignore_list.end(),
                 mol_ignore_list.begin(), mol_ignore_list.end());
        }
        else if(name == "Others" || name == "Inter")
        {
            throw std::runtime_error(
                "[error] ignore molecule do not support \"Others\" or \"Inter\".");
        }
    }

    if(ignore_table.contains("particles_within"))
    {
        if(ignore_molecule_flag)
        {
            std::cerr <<
                "[warning] ignore molecule \"Intra\" or \"Self\" was defined,"
                "so this ignore particle within bond will be ignored." << std::endl;
        }

        const auto particle_within =
            toml::find<std::map<std::string, std::size_t>>(ignore_table, "particles_within");
        for(const auto& connection : particle_within)
        {
            std::vector<std::pair<std::size_t, std::size_t>> additional_list
                = topology.ignore_list_within_edge(connection.second, connection.first);
            ignore_list.insert(ignore_list.end(),
                    additional_list.begin(), additional_list.end());
        }
    }

    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    return ignore_list;
}

const std::vector<std::pair<std::string, std::string>>
read_ignore_group(const toml::value& ignore_table)
{
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;

    if(ignore_table.contains("group"))
    {
        const auto& group = toml::find(ignore_table, "group");
        if(group.contains("inter"))
        {
            const auto& group_pairs = toml::find<toml::array>(group, "inter");
            for(const auto& group_pair : group_pairs)
            {

                const auto pair = toml::get<std::array<std::string, 2>>(group_pair);
                std::cerr << "        interaction between " << pair[0] << " and " << pair[1]
                           << " will be ignored" << std::endl;
                ignore_group_pairs.push_back({ pair[0], pair[1] });
            }
        }

        if(group.contains("intra"))
        {
            const auto&  groups = toml::find<toml::array>(group, "intra");
            for(const auto& group : groups)
            {
                const std::string group_str = toml::get<std::string>(group);
                std::cerr << "        ignore group intra in "
                          << group_str << " specified" << std::endl;

                ignore_group_pairs.push_back(std::make_pair(group_str, group_str));
            }
        }
    }

    return ignore_group_pairs;
}

const ExcludedVolumeForceFieldGenerator
read_toml_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const double eps =
        toml::find<double>(global_ff_data, "epsilon") * OpenMM::KJPerKcal; // KJPermol
    const double cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index  = Utility::find_parameter<std::size_t>(param, env, "index");
        const double      radius =
            Utility::find_parameter<double>(param, env, "radius") * OpenMM::NmPerAngstrom; // nm
        radius_vec[index] = radius;
    }

    std::cerr << "    Global        : ExcludedVolume (" << params.size()
               << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return ExcludedVolumeForceFieldGenerator(eps, cutoff, radius_vec, ignore_list,
                                             ignore_group_pairs, group_vec);
}

const WeeksChandlerAndersenForceFieldGenerator
read_toml_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec)
{
    using index_pairs_type = WeeksChandlerAndersenForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const double      sigma =
            Utility::find_parameter_either<double>(param, env, "sigma", "σ") *
            OpenMM::NmPerAngstrom; // nm
        const double      eps   =
            Utility::find_parameter_either<double>(param, env, "epsilon", "ε") *
            OpenMM::KJPerKcal; // KJ/mol
        sigma_vec[index] = sigma;
        eps_vec  [index] = eps;
    }

    std::cerr << "    Global        : WeeksChandlerAndersen ("
              << params.size() << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return WeeksChandlerAndersenForceFieldGenerator(sigma_vec, eps_vec, ignore_list,
                                                    ignore_group_pairs, group_vec);
}

const DebyeHuckelForceFieldGenerator
read_toml_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec)
{
    using index_pairs_type = DebyeHuckelForceFieldGenerator::index_pairs_type;

    const double cutoff = toml::find_or(global_ff_data, "cutoff", 5.5);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> charge_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index  = Utility::find_parameter<std::size_t>(param, env, "index");
        const double      charge = Utility::find_parameter<double>(param, env, "charge");
        charge_vec[index] = charge;
    }

    std::cerr << "    Global        : DebyeHuckel (" << params.size()
              << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return DebyeHuckelForceFieldGenerator(
            ionic_strength, temperature, cutoff, charge_vec, ignore_list,
            ignore_group_pairs, group_vec);
}

const iSoLFAttractiveForceFieldGenerator
read_toml_isolf_attractive_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec)
{
    using index_pairs_type = iSoLFAttractiveForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    std::vector<std::optional<double>> omega_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index = Utility::find_parameter<std::size_t>(param, env, "index");
        const double      sigma =
            Utility::find_parameter<double>(param, env, "sigma") * OpenMM::NmPerAngstrom; // nm
        const double      eps   =
            Utility::find_parameter<double>(param, env, "epsilon") * OpenMM::KJPerKcal; // KJPermol
        const double      omega =
            Utility::find_parameter<double>(param, env, "omega") * OpenMM::NmPerAngstrom; // nm

        sigma_vec[index] = sigma;
        eps_vec  [index] = eps;
        omega_vec[index] = omega;
    }

    std::cerr << "    Global        : iSoLFAttractive (" << params.size()
              << " found)" << std::endl;

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return iSoLFAttractiveForceFieldGenerator(sigma_vec, eps_vec, omega_vec, ignore_list,
            ignore_group_pairs, group_vec);
}

#endif // OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
