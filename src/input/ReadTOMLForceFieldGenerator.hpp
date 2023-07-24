#ifndef OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
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
#include "src/forcefield/GaussianCosineDihedralForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BaseStackingForceFieldGenerator.hpp"
#include "src/forcefield/ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2ExcludedVolumeForceFieldGenerator.hpp"
#include "src/forcefield/WeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/DebyeHuckelForceFieldGenerator.hpp"
#include "src/forcefield/iSoLFAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/UniformLennardJonesAttractiveForceFieldGenerator.hpp"
#include "src/forcefield/UniformWeeksChandlerAndersenForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2BasePairForceFieldGenerator.hpp"
#include "src/forcefield/ThreeSPN2CrossStackingForceFieldGenerator.hpp"
#include "src/forcefield/HarmonicCoMPullingForceFieldGenerator.hpp"

// -----------------------------------------------------------------------------
// read local force field

HarmonicBondForceFieldGenerator
read_toml_harmonic_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        indices.first  += offset;
        indices.second += offset;

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
        const bool use_periodic, const std::size_t ffgen_id)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;
    std::vector<double>                              sigmas;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        indices.first  += offset;
        indices.second += offset;

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
    return GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas, use_periodic, ffgen_id);
}

GoContactForceFieldGenerator
read_toml_go_contact_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{
    // TODO: enable to optimization based on cutoff
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        indices.first  += offset;
        indices.second += offset;

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
    return GoContactForceFieldGenerator(indices_vec, ks, r0s, use_periodic, ffgen_id);
}

const ThreeSPN2BondForceFieldGenerator
read_toml_3spn2_bond_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              k2s;
    std::vector<double>                              k4s;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::pair<std::size_t, std::size_t>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        indices.first  += offset;
        indices.second += offset;

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
    return ThreeSPN2BondForceFieldGenerator(indices_vec, k2s, k4s, v0s, use_periodic, ffgen_id);
}

const HarmonicAngleForceFieldGenerator
read_toml_harmonic_angle_ff_generator(
        const toml::value& local_ff_data, Topology& topology, const bool use_periodic)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 3>> indices_vec;
    std::vector<double>                     v0s;
    std::vector<double>                     ks;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 3>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        for(auto& idx : indices) { idx += offset; }

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
        const bool use_periodic, const std::size_t ffgen_id)
{
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
                Utility::find_parameter<std::array<std::size_t, 3>>(param, env, "indices");
            const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
            for(auto& idx : indices) { idx += offset; }

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
               indices_vec, ks, spline_table, aa_type, use_periodic, ffgen_id);
}

GaussianDihedralForceFieldGenerator
read_toml_gaussian_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     sigmas;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        for(auto& idx : indices) { idx += offset; }

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
            indices_vec, ks, theta0s, sigmas, use_periodic, ffgen_id);
}

const CosineDihedralForceFieldGenerator
read_toml_cosine_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     ns;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        for(auto& idx : indices) { idx += offset; }

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
            indices_vec, ks, theta0s, ns, use_periodic, ffgen_id);
}

const GaussianCosineDihedralForceFieldGenerator
read_toml_gaussian_cosine_dihedral_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{

    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks_gauss;
    std::vector<double>                     ks_cos;
    std::vector<double>                     theta0s_gauss;
    std::vector<double>                     theta0s_cos;
    std::vector<double>                     sigmas;
    std::vector<double>                     ns;

    for(const auto& param : params)
    {
        auto indices =
            Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
        const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        for(auto& idx : indices) { idx += offset; }

        for(auto idx : indices)
        {
            if(topology.size() <= idx)
            {
                throw std::runtime_error("[error] read_toml_gaussian_cosine_dihedral_ff_generator : index "+std::to_string(idx)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
            }
        }
        const auto& pot1 = Utility::find_parameter<toml::value>(param, env, "Gaussian");
        const auto& pot2 = Utility::find_parameter<toml::value>(param, env, "Cosine");

        const double k_gauss  =
            Utility::find_parameter<double>(pot1, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0_gauss =
            Utility::find_parameter<double>(pot1, env, "v0"); // radian
        const double sigma =
            Utility::find_parameter_either<double>(pot1, env, "sigma", "σ"); // radiuns

        const double k_cos  =
            Utility::find_parameter<double>(pot2, env, "k") * OpenMM::KJPerKcal; // KJ/mol
        const double theta0_cos =
            Utility::find_parameter<double>(pot2, env, "v0") + Constant::pi;    // radian
        const double n =
            Utility::find_parameter<size_t>(pot2, env, "n");

        // Toml input file assumes that a formula of the cosine dihedral potential is
        // "k*(1 + cos(n0 * (theta - t0)", while that of 3SPNC2 DNA model is
        // "k*(1 - cos(n0 * (theta - t0)". To adjust for this, we shift theta0 by π.

        indices_vec  .push_back(indices);
        ks_gauss     .push_back(k_gauss);
        ks_cos       .push_back(k_cos);
        theta0s_gauss.push_back(theta0_gauss);
        theta0s_cos  .push_back(theta0_cos);
        sigmas       .push_back(sigma);
        ns           .push_back(n);
    }

    if(local_ff_data.contains("topology"))
    {
        topology.add_edges(indices_vec, toml::find<std::string>(local_ff_data, "topology"));
    }

    std::cerr << "    DihedralAngle : Gaussian+Cosine (" << indices_vec.size() << " found)" << std::endl;
    return GaussianCosineDihedralForceFieldGenerator(
        indices_vec, ks_gauss, ks_cos, theta0s_gauss, theta0s_cos, sigmas, ns,
        use_periodic, ffgen_id);
}

const FlexibleLocalDihedralForceFieldGenerator
read_toml_flexible_local_dihedral_ff_generator(
        const toml::value& local_ff_data, const std::pair<std::string, std::string> aa_pair_type,
        const std::array<double, 7> fourier_table, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
{
    const auto& params = toml::find<toml::array>(local_ff_data, "parameters");
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;

    if(aa_pair_type.second == "GLY") // R1-GLY case
    {
        for(const auto& param : params)
        {
            auto indices = toml::find<std::array<std::size_t, 4>>(param, "indices");
            const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
            for(auto& idx : indices) { idx += offset; }

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
                    Utility::find_parameter<std::array<std::size_t, 4>>(param, env, "indices");
                const auto offset = Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
                for(auto& idx : indices) { idx += offset; }

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
               indices_vec, ks, fourier_table, aa_pair_name, use_periodic, ffgen_id);
}

const ThreeSPN2BaseStackingForceFieldGenerator
read_toml_3spn2_base_stacking_ff_generator(
        const toml::value& local_ff_data, Topology& topology,
        const bool use_periodic, const std::size_t ffgen_id)
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
    const auto& env = local_ff_data.contains("env") ? local_ff_data.at("env") : toml::value{};

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
        indices_vec, eps_vec, r0_BS_vec, theta0_BS_vec, alpha, K_BS, use_periodic, ffgen_id);
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
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic, const std::size_t ffgen_id)
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
        const std::size_t index  = Utility::find_parameter<std::size_t>(param, env, "index") +
                                   Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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
            eps, cutoff, radius_vec, ignore_list, use_periodic, ffgen_id,
            ignore_group_pairs, group_vec);
}

const ThreeSPN2ExcludedVolumeForceFieldGenerator
read_toml_3spn2_excluded_volume_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const bool use_periodic, const std::size_t ffgen_id)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

    const auto&  params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto&  env    = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const auto   eps    = ThreeSPN2ExcludedVolumePotentialParameter::epsilon;
    const auto   cutoff = toml::find_or(global_ff_data, "cutoff", 2.0);

    // Parse parameters
    //  parameters = [ # {{{
    //   {index =   0, kind = "S"},
    //   {index =   1, kind = "A"},

    std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const auto index  = Utility::find_parameter   <std::size_t>(param, env, "index") +
                            Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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
            eps, cutoff, radius_vec, ignore_list, use_periodic, ffgen_id);
}

const WeeksChandlerAndersenForceFieldGenerator
read_toml_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic, const std::size_t ffgen_id)
{
    using index_pairs_type = WeeksChandlerAndersenForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index =
            Utility::find_parameter<std::size_t>(param, env, "index") +
            Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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
            sigma_vec, eps_vec, ignore_list, use_periodic, ffgen_id,
            ignore_group_pairs, group_vec);
}

const UniformWeeksChandlerAndersenForceFieldGenerator
read_toml_uniform_weeks_chandler_andersen_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic, const std::size_t ffgen_id)
{
    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::size_t> former_participants;
    std::vector<std::size_t> latter_participants;

    for(const auto& param : params)
    {
        const std::string& particle_name = toml::find<std::string>(param, "name");
        const std::size_t  index = Utility::find_parameter<std::size_t>(param, env, "index") +
                                   Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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

    std::cerr << "    Global        : UniformWeeksChandlerAndersen - "
              << name_pair.first << "-" << name_pair.second
              <<  " (" << former_participants.size() << "-" << latter_participants.size()
              << " found)" << std::endl;

    // ignore list generation
    using index_pairs_type = UniformWeeksChandlerAndersenForceFieldGenerator::index_pairs_type;
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return UniformWeeksChandlerAndersenForceFieldGenerator(
            system_size, epsilon, sigma, former_participants, latter_participants,
            ignore_list, use_periodic, ffgen_id, ignore_group_pairs, group_vec);
}

const DebyeHuckelForceFieldGenerator
read_toml_debye_huckel_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size,
    const double ionic_strength, const double temperature, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic, const std::size_t ffgen_id)
{
    using index_pairs_type = DebyeHuckelForceFieldGenerator::index_pairs_type;

    const double cutoff = toml::find_or(global_ff_data, "cutoff", 5.5);

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> charge_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index  = Utility::find_parameter<std::size_t>(param, env, "index") +
                                   Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
        if(topology.size() <= index)
        {
            throw std::runtime_error("[error] read_toml_debye_huckel_ff_generator : index "+std::to_string(index)+" exceeds the system's largest index "+std::to_string(topology.size()-1)+".");
        }

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
            ionic_strength, temperature, cutoff, charge_vec, ignore_list, use_periodic,
            ffgen_id, ignore_group_pairs, group_vec);
}

const iSoLFAttractiveForceFieldGenerator
read_toml_isolf_attractive_ff_generator(
    const toml::value& global_ff_data, const std::size_t system_size, const Topology& topology,
    const std::vector<std::optional<std::string>>& group_vec,
    const bool use_periodic, const std::size_t ffgen_id)
{
    using index_pairs_type = iSoLFAttractiveForceFieldGenerator::index_pairs_type;

    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    std::vector<std::optional<double>> sigma_vec(system_size, std::nullopt);
    std::vector<std::optional<double>> eps_vec  (system_size, std::nullopt);
    std::vector<std::optional<double>> omega_vec(system_size, std::nullopt);
    for(const auto& param : params)
    {
        const std::size_t index = Utility::find_parameter<std::size_t>(param, env, "index") +
                                  Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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
            use_periodic, ffgen_id, ignore_group_pairs, group_vec);
}

const UniformLennardJonesAttractiveForceFieldGenerator
read_toml_uniform_lennard_jones_attractive_ff_generator(
        const toml::value& global_ff_data, const std::size_t system_size,
        const double sigma, const double epsilon,
        const std::pair<std::string, std::string>& name_pair, const Topology& topology,
        const std::vector<std::optional<std::string>>& group_vec,
        const bool use_periodic, const std::size_t ffgen_id)
{
    const auto& params = toml::find<toml::array>(global_ff_data, "parameters");
    const auto& env = global_ff_data.contains("env") ? global_ff_data.at("env") : toml::value{};

    const double cutoff = Utility::find_parameter_or<double>(global_ff_data, env, "cutoff", 2.5);

    std::vector<std::size_t> former_participants;
    std::vector<std::size_t> latter_participants;

    for(const auto& param : params)
    {
        const std::string& particle_name = toml::find<std::string>(param, "name");
        const std::size_t  index = Utility::find_parameter<std::size_t>(param, env, "index") +
                                   Utility::find_parameter_or<std::size_t>(param, env, "offset", 0);
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

    std::cerr << "    Global        : UniformLennardJonesAttractive - "
              << name_pair.first << "-" << name_pair.second
              << " (" << former_participants.size() << "-" << latter_participants.size()
              << " found)" << std::endl;

    // ignore list generation
    using index_pairs_type = UniformLennardJonesAttractiveForceFieldGenerator::index_pairs_type;
    index_pairs_type ignore_list;
    std::vector<std::pair<std::string, std::string>> ignore_group_pairs;
    if(global_ff_data.contains("ignore"))
    {
        const auto& ignore = toml::find(global_ff_data, "ignore");
        ignore_list        = read_ignore_molecule_and_particles_within(ignore, topology);
        ignore_group_pairs = read_ignore_group(ignore);
    }

    return UniformLennardJonesAttractiveForceFieldGenerator(
        system_size, epsilon, sigma, cutoff, former_participants, latter_participants,
        ignore_list, use_periodic, ffgen_id, ignore_group_pairs, group_vec);
}

template<typename PotentialParameterType>
const ThreeSPN2BasePairForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_base_pair_ff_generator(
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> base_pair,
        const bool use_periodic, const std::size_t ffgen_id
)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

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
        use_periodic, ffgen_id);
}

template<typename PotentialParameterType>
const ThreeSPN2CrossStackingForceFieldGenerator<PotentialParameterType>
read_toml_3spn2_cross_stacking_ff_generator(
        const toml::value& global_ff_data, Topology& topology,
        const std::pair<std::string, std::string> bp_kind, const std::string& strand_kind,
        const bool use_periodic, const std::size_t ffgen_id
)
{
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

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
        use_periodic, ffgen_id);
}

// -----------------------------------------------------------------------------
// read external force field

HarmonicCoMPullingForceFieldGenerator
read_toml_harmonic_com_pulling_ff_generator(
        const toml::value& external_ff_param, const bool use_periodic,
        const toml::value& env, const std::size_t ffgen_id)
{
    std::vector<int> first_group_indices;
    std::vector<int> second_group_indices;

    const double k  = Utility::find_parameter<double>(external_ff_param, env, "k");
    const double v0 = Utility::find_parameter<double>(external_ff_param, env, "v0");

    const auto& indices_pair = toml::find<toml::array>(external_ff_param, "indices_pair");
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
            k, v0, first_group_indices, second_group_indices, use_periodic, ffgen_id);
}

#endif // OPEN_AICG2_PLUS_READ_TOML_FORCE_FIELD_GENERATOR_HPP
