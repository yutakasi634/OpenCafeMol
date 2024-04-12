#ifndef OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include "Utility.hpp"
#include "src/util/Utility.hpp"
#include "src/util/Constants.hpp"
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
#include "src/forcefield/ProteinDNANonSpecificForceFieldGenerator.hpp"
#include "src/forcefield/PullingForceFieldGenerator.hpp"
#include "src/forcefield/PositionRestraintForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicCoMPullingForceFieldGenerator.hpp"

// -----------------------------------------------------------------------------
// read local force field

HarmonicBondForceFieldGenerator
read_toml_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        if(topology.size() <= indices.first)
        {
            throw std::runtime_error("[error] read_toml_harmonic_bond_ff_generator : index "+std::to_string(indices.first)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }
        else if(topology.size() <= indices.second)
        {
            throw std::runtime_error("[error] read_toml_harmonic_bond_ff_generator : index "+std::to_string(indices.second)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double v0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // nm
        // Toml input file assume the potential formula of HarmonicBond is
        // "k*(r - r0)^2", but OpenMM HarmonicBond is "1/2*k*(r - r0)^2".
        // So we needs double the interaction coefficient `k`.
        const double k =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm * 2.0; // KJ/(mol nm^2)

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : Harmonic (" << indices_vec.size() << " found)" << std::endl;
    return HarmonicBondForceFieldGenerator(indices_vec, v0s, ks, use_periodic);
}

GaussianBondForceFieldGenerator
read_toml_gaussian_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;
    std::vector<double>                              sigmas;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        if(topology.size() <= indices.first)
        {
            throw std::runtime_error("[error] read_toml_gaussian_bond_ff_generator : index "+std::to_string(indices.first)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }
        else if(topology.size() <= indices.second)
        {
            throw std::runtime_error("[error] read_toml_gaussian_bond_ff_generator : index "+std::to_string(indices.second)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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
    return GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas, use_periodic);
}

GoContactForceFieldGenerator
read_toml_go_contact_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    // TODO: enable to optimization based on cutoff
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        if(topology.size() <= indices.first)
        {
            throw std::runtime_error("[error] read_toml_go_contact_ff_generator : index "+std::to_string(indices.first)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }
        else if(topology.size() <= indices.second)
        {
            throw std::runtime_error("[error] read_toml_go_contact_ff_generator : index "+std::to_string(indices.second)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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
    return GoContactForceFieldGenerator(indices_vec, ks, r0s, use_periodic);
}

// TODO: enable to use offset
const ThreeSPN2BondForceFieldGenerator
read_toml_3spn2_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              k2s;
    std::vector<double>                              k4s;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        const double v0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // nm
        // The potential energy of the 3SPN2 bond is k (r-v0)^2 + 100*k*(r-v0)^4,
        // where the parameter, k ,is shared in the harmonic and quartic terms.
        // Unit of k for harmonic and quartic terms are kJ/(mol nm^2) and kJ/(mol nm^4),
        // respectively. Thus, the coefficient, k, for the quartic term is 100 times larger than
        // the coefficient for the harmonic term.
        const double k2 =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;  // KJ/(mol nm^2)

        const double k4 =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;  // KJ/(mol nm^4)

        indices_vec.push_back(indices);
        k2s        .push_back(k2);
        k4s        .push_back(k4);
        v0s        .push_back(v0);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    BondLength    : 3SPN2 (" << indices_vec.size() << " found)" << std::endl;
    return ThreeSPN2BondForceFieldGenerator(indices_vec, k2s, k4s, v0s, use_periodic);
}

const HarmonicAngleForceFieldGenerator
read_toml_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     v0s;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 3>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_toml_harmonic_angle_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }

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

        // Toml input file assume the potential formula of HarmonicAngle is
        // "k*(theta - theta0)^2", but OpenMM HarmonicAngle is
        // "1/2*k*(theta - theta0)^2". So we needs double the interaction
        // coefficient `k`.
        const double k =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal * 2.0; // KJ/mol

        indices_vec.push_back(indices);
        v0s        .push_back(v0);
        ks         .push_back(k);
    }

    std::cerr << "    BondAngle     : Harmonic (" << indices_vec.size() << " found)" << std::endl;

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    return HarmonicAngleForceFieldGenerator(indices_vec, v0s, ks, use_periodic);
}

FlexibleLocalAngleForceFieldGenerator
read_toml_flexible_local_angle_ff_generator(
        const toml::value& local_ff_data, const std::string& aa_type,
        const std::array<double, 10> spline_table, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        const std::string y = toml::find<std::string>(param, "y");
        if(y.substr(3, 3) == aa_type) // y is like "y1_PHE"
        {
            auto indices =
                Utility::find_parameter<std::array<std::size_t, 3>>(
                        param, env, "indices");
            const auto offset =
                Utility::find_parameter_or<toml::value>(
                        param, env, "offset", toml::value(0));
            add_offset(indices, offset);

            for(auto idx : indices)
            {
                if(topology.size() <= idx)
                {
                    throw std::runtime_error("[error] read_toml_flexible_local_angle_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
                }
            }

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
               indices_vec, ks, spline_table, aa_type, use_periodic);
}

GaussianDihedralForceFieldGenerator
read_toml_gaussian_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env
        = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     sigmas;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_toml_gaussian_dihedral_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }

        const double k  =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0 =
            Utility::find_parameter<double>(param, env, "v0"); // radian
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
    return GaussianDihedralForceFieldGenerator(
            indices_vec, ks, theta0s, sigmas, use_periodic);
}

const CosineDihedralForceFieldGenerator
read_toml_cosine_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     ns;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(
                    param, env, "indices");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(indices, offset);

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_toml_cosine_dihedral_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }

        const double k  =
            Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0 =
            Utility::find_parameter<double>(param, env, "v0") + Constant::pi;    // radian
        const double n =
            Utility::find_parameter<size_t>(param, env, "n");

        // Toml input file assumes that a formula of the cosine dihedral potential is
        // "k*(1 + cos(n0 * (theta - t0)", while that of 3SPNC2 DNA model is
        // "k*(1 - cos(n0 * (theta - t0)". To adjust for this, we shift theta0 by π.

        indices_vec.push_back(indices);
        ks         .push_back(k);
        theta0s    .push_back(theta0);
        ns         .push_back(n);


    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    DihedralAngle : Cosine (" << indices_vec.size() << " found)" << std::endl;
    return CosineDihedralForceFieldGenerator(
            indices_vec, ks, theta0s, ns, use_periodic);
}

const FlexibleLocalDihedralForceFieldGenerator
read_toml_flexible_local_dihedral_ff_generator(
        const toml::value& local_ff_data, const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table, Topology& topology,
        const bool use_periodic)
{
    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;

    if(aa_pair_type.second == "GLY") // R1-GLY case
    {
        for(const auto& param : params)
        {
            auto indices =
                toml::find<std::array<std::size_t, 4>>(param, "indices");
            const auto offset =
                Utility::find_parameter_or<toml::value>(
                        param, env, "offset", toml::value(0));
            add_offset(indices, offset);

            const std::string coef = toml::find<std::string>(param, "coef");
            const double      k    =
                Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal;

            if(coef.substr(4, 3) == "GLY")
            {
                indices_vec.push_back(indices);
                ks         .push_back(k);
            }
        }
    }
    else if(aa_pair_type.second == "PRO")
    {
        if(aa_pair_type.first == "GLY") // GLY-PRO case
        {
            for(const auto& param : params)
            {
                auto indices =
                    Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
                const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
                for(auto& idx : indices) { idx += offset; }

                const std::string coef    = toml::find<std::string>(param, "coef");
                const double      k       =
                    Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal;

                if(coef == "GLY-PRO")
                {
                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                }
            }
        }
        else
        {
            for(const auto& param : params)
            {
                auto indices =
                    Utility::find_parameter<std::array<std::size_t, 4>>(
                            param, env, "indices");
                const auto offset =
                    Utility::find_parameter_or<toml::value>(
                            param, env, "offset", toml::value(0));
                add_offset(indices, offset);

                const std::string coef    = toml::find<std::string>(param, "coef");
                const double      k       =
                    Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal;

                if(coef.substr(4, 3) == "PRO") // R2-PRO case
                {
                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                }
            }
        }
    }
    else // R1-R3 case
    {
        for(const auto& param : params)
        {
            auto indices =
                Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
            const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
            for(auto& idx : indices) { idx += offset; }

            const std::string coef    = toml::find<std::string>(param, "coef");
            const double      k       =
                Utility::find_parameter<double>(param, env, "k") * OpenMM::KJPerKcal;

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
    }

    for(auto& indices : indices_vec)
    {
        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_toml_flexible_local_dihedral_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    const std::string aa_pair_name =
        aa_pair_type.first + "-" + aa_pair_type.second;
    std::cerr << "    DihedralAngle : FlexibleLocalDihedral - "
              << aa_pair_name << " (" << indices_vec.size() << " found)" << std::endl;
    return FlexibleLocalDihedralForceFieldGenerator(
               indices_vec, ks, fourier_table, aa_pair_name, use_periodic);
}

// TODO: enable to use offset
const ThreeSPN2BaseStackingForceFieldGenerator
read_toml_3spn2_base_stacking_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic)
{
    // [[forcefields.local]]
    // interaction = "3SPN2BaseStacking"
    // potential   = "3SPN2C"
    // topology    = "nucleotide"
    // parameters = [
    //     {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    check_keys_available(local_ff_data,
            {"interaction", "potential", "topology", "parameters", "env"});

    const auto pot = toml::find<std::string>(local_ff_data, "potential");
    if (! (pot == "3SPN2" || pot == "3SPN2C"))
    {
        throw std::runtime_error(
            "[error] invalid potential " + pot + " fond."
            "Expected value is one of the following."
            "- \"3SPN2\" : The general 3SPN2 parameter set."
            "- \"3SPN2C\": The parameter set optimized to reproduce sequence-dependent curveture of dsDNA."
        );
    }

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env =
        local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

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
    //     {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
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

    // Parameter assignment
    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double> eps_vec;
    std::vector<double> r0_BS_vec;
    std::vector<double> theta0_BS_vec;
    double alpha = 1.0;
    double K_BS  = 1.0;

    for (std::size_t i=1; i<nucleotide_idxs.size(); ++i)
    {
        const auto& Base5 = nucleotide_idxs.at(i-1);
        const auto& Base3 = nucleotide_idxs.at(i);
        const std::string base_kind = Base5.base + Base3.base;

        if(Base5.strand != Base3.strand)
        {
            continue; // if the strands are different, there is no base stacking
        }
        assert(Base3.base != "X");
        assert(Base5.base != "X");

        const std::array<std::size_t, 3> indices{{Base5.S, Base5.B, Base3.B}};

        if (pot == "3SPN2")
        {
            using parameter_type = ThreeSPN2BaseStackingPotentialParameter;
            const double eps       = parameter_type::epsilon_BS.at(base_kind); // kJ/mol
            const double r0_BS     = parameter_type::r0_BS     .at(base_kind) * OpenMM::NmPerAngstrom; // nm
            const double theta0_BS = parameter_type::theta0_BS .at(base_kind) * OpenMM::RadiansPerDegree; // radian
            alpha                  = parameter_type::alpha_BS / OpenMM::NmPerAngstrom; // nm^{-1}
            K_BS                   = parameter_type::K_BS;

            indices_vec  .push_back(indices);
            eps_vec      .push_back(eps);
            r0_BS_vec    .push_back(r0_BS);
            theta0_BS_vec.push_back(theta0_BS);
        }
        else if (pot == "3SPN2C")
        {
            using parameter_type = ThreeSPN2CBaseStackingPotentialParameter;
            const double eps       = parameter_type::epsilon_BS.at(base_kind); // kJ/mol
            const double r0_BS     = parameter_type::r0_BS     .at(base_kind) * OpenMM::NmPerAngstrom;    // nm
            const double theta0_BS = parameter_type::theta0_BS .at(base_kind) * OpenMM::RadiansPerDegree; // radian
            alpha                  = parameter_type::alpha_BS / OpenMM::NmPerAngstrom; // nm^{-1}
            K_BS                   = parameter_type::K_BS;

            indices_vec  .push_back(indices);
            eps_vec      .push_back(eps);
            r0_BS_vec    .push_back(r0_BS);
            theta0_BS_vec.push_back(theta0_BS);
        }
    }

    // Creating topology
    if (local_ff_data.contains("topology"))
    {
        for (std::size_t i=1; i<nucleotide_idxs.size(); ++i)
        {
            using edge_type = std::pair<std::size_t, std::size_t>;
            std::vector<edge_type> edges;

            const auto& Base5 = nucleotide_idxs.at(i-1);
            const auto& Base3 = nucleotide_idxs.at(i);
            constexpr auto nil = Nucleotide::nil();

            if (Base5.strand != Base3.strand) {continue;}

            edges.push_back(std::make_pair(Base5.S, Base3.S));
            edges.push_back(std::make_pair(Base5.S, Base3.B));
            edges.push_back(std::make_pair(Base5.B, Base3.S));
            edges.push_back(std::make_pair(Base5.B, Base3.B));

            if (Base5.P != nil)
            {
                edges.push_back(std::make_pair(Base5.P, Base3.S));
                edges.push_back(std::make_pair(Base5.P, Base3.B));
            }
            if (Base3.P != nil)
            {
                edges.push_back(std::make_pair(Base5.S, Base3.P));
                edges.push_back(std::make_pair(Base5.B, Base3.P));
            }
            if (Base5.P != nil && Base3.P != nil)
            {
                edges.push_back(std::make_pair(Base5.P, Base3.P));
            }
            topology.add_edges(edges, toml::find<std::string>(local_ff_data, "topology"));
        }
    }

    if (pot == "3SPN2")
    {
        std::cerr << "    3SPN2BaseStacking    : 3SPN2 ("  << indices_vec.size() << " found)" << std::endl;
    }
    else if (pot == "3SPN2C")
    {
        std::cerr << "    3SPN2BaseStacking    : 3SPN2C (" << indices_vec.size() << " found)" << std::endl;
    }

    return ThreeSPN2BaseStackingForceFieldGenerator(
        indices_vec, eps_vec, r0_BS_vec, theta0_BS_vec, alpha, K_BS, use_periodic);
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
            std::cerr << "\033[33m[warning]\033[m"
                      << " ignore molecule \"Intra\" or \"Self\" was defined,"
                         "so this ignore particle within bond will be ignored."
                      << std::endl;
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
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "env",
             "cutoff", "parameters", "epsilon"});

    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const double eps =
        toml::find<double>(global_ff_data, "epsilon") * OpenMM::KJPerKcal; // KJPermol
    const double cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_excluded_volume_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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

    return ExcludedVolumeForceFieldGenerator(
            eps, cutoff, radius_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

const ThreeSPN2ExcludedVolumeForceFieldGenerator
read_toml_3spn2_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const Topology& topology, const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "cutoff", "parameters", "env"});

    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const auto&  params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto&  env
        = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const auto   eps    = ThreeSPN2ExcludedVolumePotentialParameter::epsilon;
    const auto   cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    // Parse parameters
    //  parameters = [ # {{{
    //   {index =   0, kind = "S"},
    //   {index =   1, kind = "A"},

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index  =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        const auto kind   = toml::find<std::string>(param, "kind");
        radius_vec[index] = ThreeSPN2ExcludedVolumePotentialParameter::sigma.at(kind);
    }

    std::cerr << "    Global        : ExcludedVolume 3SPN2 (" << params.size()
               << " found)" << std::endl;

    // ignore list generation for bonded interactions and neighboring nucleotides
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    // ignore list generation for base pairing interaction
    for(size_t idx=0; idx<params.size()-1; ++idx)
    {
        for(size_t jdx=idx+1; jdx<params.size(); ++jdx)
        {
            const auto kind_i = toml::find<std::string>(params[idx], "kind");
            const auto kind_j = toml::find<std::string>(params[jdx], "kind");

            // complement
            if ((kind_i == "A" && kind_j == "T") || (kind_i == "T" && kind_j == "A") ||
                (kind_i == "G" && kind_j == "C") || (kind_i == "C" && kind_j == "G"))
            {
                const auto index_i =
                    Utility::find_parameter   <std::size_t>(params[idx], env, "index") +
                    Utility::find_parameter_or<std::size_t>(params[idx], env, "offset", 0);
                const auto index_j =
                    Utility::find_parameter   <std::size_t>(params[jdx], env, "index") +
                    Utility::find_parameter_or<std::size_t>(params[jdx], env, "offset", 0);

                if (index_i < index_j)
                {
                    ignore_list.push_back(std::make_pair(index_i, index_j));
                }
                else
                {
                    ignore_list.push_back(std::make_pair(index_j, index_i));
                }
            }
        }
    }

    // remove duplication
    std::sort(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(std::unique(ignore_list.begin(), ignore_list.end()), ignore_list.end());

    return ThreeSPN2ExcludedVolumeForceFieldGenerator(
            eps, cutoff, radius_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

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
        const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "parameters", "env"});

    using index_pairs_type = WeeksChandlerAndersenForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_weeks_chandler_andersen_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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

    return WeeksChandlerAndersenForceFieldGenerator(
            sigma_vec, eps_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

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
        const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "parameters", "table", "env"});

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::size_t> former_participants;
    std::vector<std::size_t> latter_participants;

    for(const auto& param : params)
    {
        const std::string& particle_name = toml::find<std::string>(param, "name");
        std::size_t  index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_uniform_weeks_chandler_andersen_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        if(particle_name == name_pair.first)
        {
            former_participants.push_back(index);
        }
        if(particle_name == name_pair.second)
        {
            latter_participants.push_back(index);
        }
    }

    std::cerr << "        "
              << name_pair.first << "-" << name_pair.second
              <<  " (" << former_participants.size() << "-" << latter_participants.size()
              << " found)" << std::endl;

    return UniformWeeksChandlerAndersenForceFieldGenerator(
            system_size, epsilon, sigma, former_participants, latter_participants,
            ignore_list, use_periodic, ignore_group_pairs, group_vec);
}

// DebyeHuckel input is like below
// [[forcefields.global]]
// interaction = "Pair"
// potential   = "DebyeHuckel"
// ignore.particles_within.bond = 3
// parameters = [ # {{{
//     {index =   2, charge = -0.6},
//     # ...
// ]
const DebyeHuckelForceFieldGenerator
read_toml_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "cutoff", "parameters", "env"});

    using index_pairs_type = DebyeHuckelForceFieldGenerator::index_pairs_type;

    const double cutoff = toml::find_or(global_ff_data, "cutoff", 5.5);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> charge_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index  =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_debye_huckel_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double charge = Utility::find_parameter<double>(param, env, "charge");
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
            ionic_strength, temperature, cutoff, charge_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

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
    const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "parameters", "env"});

    using index_pairs_type = iSoLFAttractiveForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    std::vector<std::optional<double>> omega_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_isolf_attractive_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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
            use_periodic, ignore_group_pairs, group_vec);
}


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
        const bool use_periodic)
{
    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "table", "cutoff", "parameters", "env"});

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const double cutoff = Utility::find_parameter_or<double>(global_ff_data, env, "cutoff", 2.5);

    std::vector<std::size_t> former_participants;
    std::vector<std::size_t> latter_participants;

    for(const auto& param : params)
    {
        const std::string& particle_name = toml::find<std::string>(param, "name");
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_uniform_lennard_jones_attractive_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        if(particle_name == name_pair.first)
        {
            former_participants.push_back(index);
        }
        if(particle_name == name_pair.second)
        {
            latter_participants.push_back(index);
        }
    }

    std::cerr << "        "
              << name_pair.first << "-" << name_pair.second
              << " (" << former_participants.size() << "-" << latter_participants.size()
              << " found)" << std::endl;

    return UniformLennardJonesAttractiveForceFieldGenerator(
        system_size, epsilon, sigma, cutoff, former_participants, latter_participants,
        ignore_list, use_periodic, ignore_group_pairs, group_vec);
}

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
    // parameters  = [ # {{{
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
    // parameters  = [ # {{{
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
        const bool use_periodic)
{
    using index_pairs_type = LennardJonesAttractiveForceFieldGenerator::index_pairs_type;

    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "cutoff", "parameters", "env"});

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env =
        global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};
    const double cutoff_ratio =
        Utility::find_parameter_or<double>(global_ff_data, env, "cutoff", 2.5);

    std::vector<std::optional<double>> epsilon_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> sigma_vec  (system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index") +
            Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_lennard_jones_attractive_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double sigma =
            Utility::find_parameter<double>(param, env, "sigma") *
            OpenMM::NmPerAngstrom; // nm
        const double epsilon =
            Utility::find_parameter<double>(param, env, "epsilon") *
            OpenMM::KJPerKcal; // KJPermol

        epsilon_vec[index] = epsilon;
        sigma_vec  [index] = sigma;
    }

    std::cerr << "    Global        : LennardJonesAttractive (" << params.size()
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

    return LennardJonesAttractiveForceFieldGenerator(cutoff_ratio,
            epsilon_vec, sigma_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}


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
        const bool use_periodic)
{
    using index_pairs_type = LennardJonesRepulsiveForceFieldGenerator::index_pairs_type;

    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "cutoff", "parameters", "env"});

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env =
        global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};
    const double cutoff_ratio =
        Utility::find_parameter_or<double>(global_ff_data, env, "cutoff", 2.5);

    std::vector<std::optional<double>> epsilon_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> sigma_vec  (system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index") +
            Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);

        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_lennard_jones_repulsive_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

        const double sigma =
            Utility::find_parameter<double>(param, env, "sigma") *
            OpenMM::NmPerAngstrom; // nm
        const double epsilon =
            Utility::find_parameter<double>(param, env, "epsilon") *
            OpenMM::KJPerKcal; // KJPermol

        epsilon_vec[index] = epsilon;
        sigma_vec  [index] = sigma;
    }

    std::cerr << "    Global        : LennardJonesRepulsive (" << params.size()
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

    return LennardJonesRepulsiveForceFieldGenerator(cutoff_ratio,
            epsilon_vec, sigma_vec, ignore_list, use_periodic,
            ignore_group_pairs, group_vec);
}

CombinatorialGoContactForceFieldGenerator
read_toml_combinatorial_go_contact_ff_generator(
        const double cutoff_ratio, const toml::value& contacts_info,
        const toml::value& env, const Topology& topology,
        const std::vector<std::pair<std::size_t, std::size_t>>& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>>& ignore_group_pairs,
        const std::vector<std::optional<std::string>>& group_vec)
{
    auto indices_pair =
        Utility::find_parameter<
            std::pair<std::vector<std::size_t>,
                      std::vector<std::size_t>>>(contacts_info, env, "indices_pair");

    const auto offset =
        Utility::find_parameter_or<std::size_t>(contacts_info, env, "offset", 0);
    const double k  =
        Utility::find_parameter<double>(contacts_info, env, "k") * OpenMM::KJPerKcal; // KJ/mol
    const double r0 =
        Utility::find_parameter<double>(contacts_info, env, "v0") * OpenMM::NmPerAngstrom; // nm

    std::vector<std::size_t> first_indices, second_indices;
    for(std::size_t idx : indices_pair.first)
    {
        first_indices.push_back(idx + offset);
    }
    for(std::size_t idx : indices_pair.second)
    {
        second_indices.push_back(idx + offset);
    }

    return CombinatorialGoContactForceFieldGenerator(
            k, r0, cutoff_ratio, indices_pair.first, indices_pair.second,
            ignore_list, use_periodic, ignore_group_pairs, group_vec);
}

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
        const bool use_periodic)
{
    using index_pairs_type =
        CombinatorialGoContactForceFieldGenerator::index_pairs_type;

    check_keys_available(global_ff_data,
            {"interaction", "potential", "ignore", "cutoff", "parameters", "env"});

    const auto& parameters = toml::find<toml::array>(global_ff_data, "parameters");
    std::cerr << "    Global       : CombinatorialGoContact (" << parameters.size()
              << " found)" << std::endl;

    const auto& env =
        global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const double cutoff_ratio = toml::find_or(global_ff_data, "cutoff", 1.8);

    // ignore list generation
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list =
            read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    std::vector<CombinatorialGoContactForceFieldGenerator> ff_gens;
    for(const auto& contacts_info : parameters)
    {
        CombinatorialGoContactForceFieldGenerator ff_gen =
            read_toml_combinatorial_go_contact_ff_generator(
                cutoff_ratio, contacts_info, env, topology, ignore_list,
                use_periodic, ignore_group_pairs, group_vec);
        ff_gens.push_back(ff_gen);
    }

    return ff_gens;
}

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
    const  toml::value& global_ff_data, const bool use_periodic
)
{
    check_keys_available(global_ff_data,
    { "interaction", "potential", "env", "sigma", "delta", "cutoff", "energy_unit", "parameters"});

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const double sigma =
        Utility::find_parameter_or<double>(global_ff_data, env, "sigma", 1.0) * OpenMM::NmPerAngstrom; // nm;

    const double delta =
        Utility::find_parameter_or<double>(global_ff_data, env, "delta", 0.17453); // radian

    const double cutoff_ratio =
        Utility::find_parameter_or<double>(global_ff_data, env, "cutoff", 10.0);

    const double energy_unit  =
        Utility::find_parameter_or<double>(global_ff_data, env, "energy_unit", 0.593) * OpenMM::KJPerKcal; // kJ


    // Vector for indices
    std::vector<std::array<std::size_t, 2>> indices_dna;     // 0: phos, 1: sugar3
    std::vector<std::array<std::size_t, 3>> indices_protein; // 0: CA, 1: CA_N, 2: CA_C

    // Vector for PDNS parameters
    std::vector<double> ks;
    std::vector<double> r0s;
    std::vector<double> theta0s;
    std::vector<double> phi0s;

    for(const auto& param: params)
    {
        std::size_t index  =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);

        const auto kind = toml::find<std::string>(param, "kind");

        if (kind == "DNA")
        {
            const::size_t idx_sugar3 = Utility::find_parameter<std::size_t>(param, env, "S3");
            indices_dna.push_back({index, idx_sugar3});
        }

        if (kind == "Protein")
        {
            const::size_t idx_calpha_n = Utility::find_parameter<std::size_t>(param, env, "PN");
            const::size_t idx_calpha_c = Utility::find_parameter<std::size_t>(param, env, "PC");

            indices_protein.push_back({index, idx_calpha_n, idx_calpha_c});

            const double k =
                Utility::find_parameter<double>(param, env, "k") * energy_unit; // kJ

            const double r0 =
                Utility::find_parameter<double>(param, env, "r0") * OpenMM::NmPerAngstrom; // nm

            const double theta0 =
                Utility::find_parameter<double>(param, env, "theta0") * OpenMM::RadiansPerDegree; // radian

            const double phi0 =
                Utility::find_parameter<double>(param, env, "phi0") * OpenMM::RadiansPerDegree; // radian

            ks.push_back(k);
            r0s.push_back(r0);
            theta0s.push_back(theta0);
            phi0s.push_back(phi0);
        }
    }

    return ProteinDNANonSpecificForceFieldGenerator(
        indices_dna, indices_protein, sigma, delta, cutoff_ratio,
        ks, r0s, theta0s, phi0s, use_periodic);
}

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
        const bool use_periodic)
{
    using parameter_type = PullingForceFieldGenerator::parameter_type;

    check_keys_available(external_ff_data, {"interaction", "parameters", "env"});

    const auto& params = toml::find<toml::array>(external_ff_data, "parameters");
    const auto& env =
        external_ff_data.contains("env") ? external_ff_data.at("env") : toml::value{};

    std::vector<parameter_type> idx_force_vec;
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);
        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_pulling_ff_generator : index " +
                std::to_string(index) + " exceeds the system's largest index " +
                std::to_string(topology.size()-1) + ".");
        }
        const std::array<double, 3> force_kcal_vec =
            Utility::find_parameter<std::array<double, 3>>(param, env, "force"); // [kcal/(mol Å)]
        std::array<double, 3> force_kj_vec;
        for(std::size_t idx = 0; idx < 3; idx++)
        {
            force_kj_vec[idx] =
                force_kcal_vec[idx] * OpenMM::KJPerKcal * OpenMM::AngstromsPerNm;
        }
        idx_force_vec.push_back({index, force_kj_vec});
    }

    std::cerr << "    External      : Pulling (" << params.size()
              << " found)" << std::endl;

    return PullingForceFieldGenerator(idx_force_vec, use_periodic);
}

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
        const toml::value& external_ff_data, const Topology& topology)
{
    check_keys_available(external_ff_data,
            {"interaction", "potential", "parameters", "env"});

    const auto& params = toml::find<toml::array>(external_ff_data, "parameters");
    const auto& env =
        external_ff_data.contains("env") ? external_ff_data.at("env") : toml::value{};

    std::vector<std::size_t>           indices;
    std::vector<std::array<double, 3>> positions;
    std::vector<double>                ks;
    std::vector<double>                v0s;
    for(const auto& param : params)
    {
        std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index");
        const auto offset =
            Utility::find_parameter_or<toml::value>(
                    param, env, "offset", toml::value(0));
        add_offset(index, offset);
        if(topology.size() <= index)
        {
            throw std::runtime_error(
                "[error] read_toml_position_restraint_ff_generator : index " +
                std::to_string(index) + " exceeds the system's largest index " +
                std::to_string(topology.size()-1) + ".");
        }
        std::array<double, 3> position =
            Utility::find_parameter<std::array<double, 3>>(param, env, "position"); // [Å]
        const double k =
            Utility::find_parameter<double>(param, env, "k" ) * OpenMM::KJPerKcal *
            OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm; // [kJ/(mol nm^2)]
        const double v0 =
            Utility::find_parameter<double>(param, env, "v0") * OpenMM::NmPerAngstrom; // [Nm]
        for(std::size_t idx = 0; idx < 3; idx++)
        {
            position[idx] = position[idx] * OpenMM::NmPerAngstrom; // [Nm]
        }

        indices  .push_back(index);
        positions.push_back(position);
        ks       .push_back(k);
        v0s      .push_back(v0);
    }

    std::cerr << "    External      : PositionRestraint (" << params.size()
              << " found)" << std::endl;

    return PositionRestraintForceFieldGenerator(indices, positions, ks, v0s);
}

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
        const toml::value& env)
{
    check_keys_available(external_ff_data,
            {"interaction", "potential", "parameters", "env"});

    std::vector<int> first_group_indices;
    std::vector<int> second_group_indices;

    const double k  = Utility::find_parameter<double>(external_ff_data, env, "k");
    const double v0 = Utility::find_parameter<double>(external_ff_data, env, "v0");

    const auto& indices_pair = toml::find<toml::array>(external_ff_data, "indices_pair");
    for(const auto& idx : indices_pair.at(0).as_array())
    {
        first_group_indices.push_back(toml::get<std::size_t>(idx));
    }

    for(const auto& idx : indices_pair.at(1).as_array())
    {
        second_group_indices.push_back(toml::get<std::size_t>(idx));
    }

    std::cerr << "    External      : HarmonicCoMPulling"
              << " (" << first_group_indices.size() << "-" << second_group_indices.size()
              << " found)" << std::endl;

    return HarmonicCoMPullingForceFieldGenerator(
            k, v0, first_group_indices, second_group_indices, use_periodic);
}

#endif // OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
