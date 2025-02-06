#include "ReadGenesisForceFieldGenerator.hpp"

#include "src/util/Utility.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_bonds_section(
        const std::vector<std::string>& bonds_data, Topology& topology,
        const bool use_periodic)
{

    struct bond_params
    {
        const std::pair<std::size_t, std::size_t> indices;
        const double                              equilibrium_value;
        const double                              force_coef;

        bond_params(const std::pair<std::size_t, std::size_t> idxs,
                const double eq_val, const double f_coef)
            : indices(idxs), equilibrium_value(eq_val), force_coef(f_coef)
        {}
    };

    std::vector<bond_params> type1_bond_params; // harmonic_bond
    std::vector<bond_params> type21_bond_params;

    for(auto& bond_line : bonds_data)
    {
        const std::size_t idx_i  = std::stoi(bond_line.substr( 0, 10)) - 1;
        const std::size_t idx_j  = std::stoi(bond_line.substr(10, 10)) - 1;
        const std::size_t f_type = std::stoi(bond_line.substr(20,  5));
        const double      eq     = std::stod(bond_line.substr(25, 42));
        const double      coef   = std::stod(bond_line.substr(43, 60)) * 0.5; // KJ/(mol nm^2)
        if(f_type == 1)
        {
            type1_bond_params.push_back(
                    bond_params(std::make_pair(idx_i, idx_j), eq, coef));
        }
        else if(f_type == 21)
        {
            throw std::runtime_error("[error] bond type 21 is not yet supported.");
        }
        else
        {
            throw std::runtime_error(
                    "[error] unexpected bond type " + std::to_string(f_type) +
                    " was specified in `[ bonds ]` section. expected value is 1 or 21.");
        }
    }

    std::vector<std::unique_ptr<ForceFieldGeneratorBase>> ff_gen_ptrs{};

    if(type1_bond_params.size() != 0)
    {
        std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
        std::vector<double>                              eqs;
        std::vector<double>                              coefs;
        for(const auto& param : type1_bond_params)
        {
            indices_vec.emplace_back(param.indices);
            eqs        .emplace_back(param.equilibrium_value);
            coefs      .emplace_back(param.force_coef);
        }
        topology.add_edges(indices_vec, "bond");

        std::cerr << "    BondLength    : Harmonic (" <<
                  indices_vec.size() << " found)" << std::endl;

        ff_gen_ptrs.push_back(
                std::make_unique<HarmonicBondForceFieldGenerator>(
                    indices_vec, eqs, coefs, use_periodic));
    }

    return ff_gen_ptrs;
}

GaussianBondForceFieldGenerator
read_genesis_gaussian_bond_ff_generator(
        const std::vector<std::string>& angles_data, const bool use_periodic)
{
    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;
    std::vector<double>                              sigmas;

    for(const auto& angles_line : angles_data)
    {
        if(std::stoi(angles_line.substr(30, 5)) == 21) // 21 means AICG-type angle function
        {
            const std::size_t idx_i = std::stoi(angles_line.substr(0,  10)) - 1;
            const std::size_t idx_k = std::stoi(angles_line.substr(20, 10)) - 1;
            const double      v0    = std::stod(angles_line.substr(35, 15)); // nm
            const double      k     = -std::stod(angles_line.substr(50, 15)); // KJ/mol
            const double      sigma = std::stod(angles_line.substr(65, 15)); // nm

            indices_vec.push_back(std::make_pair(idx_i, idx_k));
            v0s        .push_back(v0);
            ks         .push_back(k);
            sigmas     .push_back(sigma);
        }
    }

    std::cerr << "    BondLength    : Gaussian (" << indices_vec.size() << " found)" << std::endl;
    return GaussianBondForceFieldGenerator(
            indices_vec, ks, v0s, sigmas, use_periodic);
}

GoContactForceFieldGenerator
read_genesis_go_contact_ff_generator(
        const std::vector<std::string>& pairs_data, Topology& topology,
        const bool use_periodic)
{
    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(auto& pairs_line : pairs_data)
    {
        const std::size_t idx_i = std::stoi(pairs_line.substr(0,  10)) - 1;
        const std::size_t idx_k = std::stoi(pairs_line.substr(10, 10)) - 1;
        const double      r0    = std::stod(pairs_line.substr(30, 15)); // nm
        const double      k     = std::stod(pairs_line.substr(45, 15)); // KJ/mol

        indices_vec.push_back(std::make_pair(idx_i, idx_k));
        ks         .push_back(k);
        r0s        .push_back(r0);
    }
    topology.add_edges(indices_vec, "contact");

    std::cerr << "    BondLength    : GoContact (" << indices_vec.size() << " found)" << std::endl;
    return GoContactForceFieldGenerator(indices_vec, ks, r0s, use_periodic);
}

FlexibleLocalAngleForceFieldGenerator
read_genesis_flexible_local_angle_ff_generator(
        const std::vector<std::string>& angles_data, const std::vector<std::string>& atoms_data,
        const std::string& aa_type, const bool use_periodic)
{
    std::vector<std::string> aa_type_vec;
    for(const auto& atoms_line : atoms_data)
    {
        const std::string atom_type = Utility::erase_space(atoms_line.substr(10, 5));
        aa_type_vec.push_back(atom_type);
    }

    std::vector<std::array<std::size_t, 3>> indices_vec;
    for(const auto& angles_line : angles_data)
    {
        if(std::stoi(angles_line.substr(30, 5)) == 22) // 22 means FLP-type angle function
        {
            const std::size_t idx_j = std::stoi(angles_line.substr(10, 10)) - 1;
            if(aa_type_vec[idx_j] == aa_type)
            {
                const std::size_t idx_i = std::stoi(angles_line.substr(0,  10)) - 1;
                const std::size_t idx_k = std::stoi(angles_line.substr(20, 10)) - 1;

                indices_vec.push_back({idx_i, idx_j, idx_k});
            }
        }
    }

    std::cerr << "    BondAngle     : FlexibleLocalAngle - "
              << aa_type << " (" << indices_vec.size() << " found)" << std::endl;

    return FlexibleLocalAngleForceFieldGenerator(
               indices_vec, std::vector<double>(indices_vec.size(), OpenMM::KJPerKcal),
               Constant::fla_spline_table.at(aa_type), aa_type, use_periodic);
}

GaussianDihedralForceFieldGenerator
read_genesis_gaussian_dihedral_ff_generator(
        const std::vector<std::string>& dihedrals_data,
        const bool use_periodic)
{
    std::vector<std::array<std::size_t, 4>> indices_vec;
    std::vector<double>                     ks;
    std::vector<double>                     theta0s;
    std::vector<double>                     sigmas;

    for(const auto& dihedrals_line : dihedrals_data)
    {
        if(std::stoi(dihedrals_line.substr(40, 5)) == 41) // 41 means safe-dihedral for AICG2+
        {
            const std::size_t  idx_i  = std::stoi(dihedrals_line.substr( 0, 10)) - 1;
            const std::size_t  idx_j  = std::stoi(dihedrals_line.substr(10, 10)) - 1;
            const std::size_t  idx_k  = std::stoi(dihedrals_line.substr(20, 10)) - 1;
            const std::size_t  idx_l  = std::stoi(dihedrals_line.substr(30, 10)) - 1;
            const double       theta0 = std::stod(dihedrals_line.substr(45, 15)) / 180.0 * Constant::pi; // radiuns
            const double       k      = -std::stod(dihedrals_line.substr(60, 15)); // KJ/mol
            const double       sigma  = std::stod(dihedrals_line.substr(75, 15)); // radiuns

            indices_vec.push_back({idx_i, idx_j, idx_k, idx_l});
            ks         .push_back(k);
            theta0s    .push_back(theta0);
            sigmas     .push_back(sigma);
        }
    }

    std::cerr << "    DihedralAngle : Gaussian (" << indices_vec.size() << " found)" << std::endl;
    return GaussianDihedralForceFieldGenerator(
            indices_vec, ks, theta0s, sigmas, use_periodic);
}

FlexibleLocalDihedralForceFieldGenerator
read_genesis_flexible_local_dihedral_ff_generator(
        const std::vector<std::string>& dihedrals_data,
        const std::vector<std::string>& atoms_data,
        const std::pair<std::string, std::string>& aa_type_pair,
        const bool use_periodic)
{
    std::vector<std::string> aa_type_vec;
    for(const auto& atoms_line : atoms_data)
    {
        aa_type_vec.push_back(Utility::erase_space(atoms_line.substr(10, 5)));
    }

    std::vector<std::array<std::size_t, 4>> indices_vec;
    for(const auto& dihedrals_line : dihedrals_data)
    {
        if(std::stoi(dihedrals_line.substr(40, 5)) == 52) // 52 means safe-dihedral for FLP
        {
            if(aa_type_pair.first == "R1")
            {
                const std::size_t idx_k = std::stoi(dihedrals_line.substr(20, 10)) - 1;
                if(aa_type_vec[idx_k] == "GLY")
                {
                    const std::size_t idx_i = std::stoi(dihedrals_line.substr( 0, 10)) - 1;
                    const std::size_t idx_j = std::stoi(dihedrals_line.substr(10, 10)) - 1;
                    const std::size_t idx_l = std::stoi(dihedrals_line.substr(30, 10)) - 1;
                    indices_vec.push_back({idx_i, idx_j, idx_k, idx_l});
                }
            }
            else if(aa_type_pair.first == "R2")
            {
                const std::size_t idx_j = std::stoi(dihedrals_line.substr(10, 10)) - 1;
                const std::size_t idx_k = std::stoi(dihedrals_line.substr(20, 10)) - 1;
                if(aa_type_vec[idx_j] != "GLY" && aa_type_vec[idx_k] == "PRO")
                {
                    const std::size_t idx_i = std::stoi(dihedrals_line.substr( 0, 10)) - 1;
                    const std::size_t idx_l = std::stoi(dihedrals_line.substr(30, 10)) - 1;
                    indices_vec.push_back({idx_i, idx_j, idx_k, idx_l});
                }
            }
            else if(aa_type_pair.first == "GLY" && aa_type_pair.second == "PRO")
            {
                const std::size_t idx_j = std::stoi(dihedrals_line.substr(10, 10)) - 1;
                const std::size_t idx_k = std::stoi(dihedrals_line.substr(20, 10)) - 1;
                if(aa_type_vec[idx_j] == "GLY" && aa_type_vec[idx_k] == "PRO")
                {
                    const std::size_t idx_i = std::stoi(dihedrals_line.substr( 0, 10)) - 1;
                    const std::size_t idx_l = std::stoi(dihedrals_line.substr(30, 10)) - 1;
                    indices_vec.push_back({idx_i, idx_j, idx_k, idx_l});
                }
            }
            else
            {
                const std::size_t idx_j = std::stoi(dihedrals_line.substr(10, 10)) - 1;
                const std::size_t idx_k = std::stoi(dihedrals_line.substr(20, 10)) - 1;
                if(aa_type_vec[idx_j] == aa_type_pair.first)
                {
                    const std::size_t idx_i = std::stoi(dihedrals_line.substr( 0, 10)) - 1;
                    const std::size_t idx_l = std::stoi(dihedrals_line.substr(30, 10)) - 1;
                    indices_vec.push_back({idx_i, idx_j, idx_k, idx_l});
                }
            }
        }
    }

    std::cerr << "    DihedralAngle : FlexibleLocalDihedral - " << aa_type_pair.first
              << "-" << aa_type_pair.second << " (" << indices_vec.size() << " found)"
              << std::endl;

    std::vector<double> ks = std::vector<double>(indices_vec.size(), OpenMM::KJPerKcal);
    return FlexibleLocalDihedralForceFieldGenerator(indices_vec, ks,
            Constant::fld_fourier_table.at(aa_type_pair),
            aa_type_pair.first + "-" + aa_type_pair.second, use_periodic);
}

ExcludedVolumeForceFieldGenerator
read_genesis_exv_ff_generator(const std::vector<std::string>& atomtypes_data,
        const std::vector<std::string>& atoms_data, Topology& topology,
        const bool use_periodic,
        const std::size_t ignore_particle_within_bond)
{
    std::map<std::string, double> name_rmin_map;
    std::string                   eps_str;
    for(auto& atomtypes_line : atomtypes_data)
    {
        const std::string name = Utility::erase_space(atomtypes_line.substr(0, 5));
        if(name_rmin_map.find(name) == name_rmin_map.end())
        {
            const double rmin = std::stod(atomtypes_line.substr(30, 10)) * 0.5;
            name_rmin_map.insert(std::make_pair(name, rmin));
        }
        else
        {
            throw std::runtime_error("[error] `" + name + "` key was duplicated in "
                                     "`[ atomtypes ]` section.");
        }

        if(eps_str.empty())
        {
            eps_str = Utility::erase_space(atomtypes_line.substr(41, 9));
        }
        else if(eps_str != Utility::erase_space(atomtypes_line.substr(41, 9)))
        {
            throw std::runtime_error("[error] ExcludedVolumeForceField only support uniform eps for all particle case. Some different value were defined in `[ atomtypes ]` table.");
        }
    }
    const double eps = std::stod(eps_str);

    std::vector<std::optional<double>> radius_vec;
    for(auto& atoms_line : atoms_data)
    {
        const std::string aa_type = Utility::erase_space(atoms_line.substr(10, 5));
        radius_vec.push_back(name_rmin_map.at(aa_type));
    }

    std::cerr << "    Global        : ExcludedVolume (" << radius_vec.size()
              << " found)"                              << std::endl;

    // generate ignore list
    ExcludedVolumeForceFieldGenerator::index_pairs_type ignore_list =
        topology.ignore_list_within_edge(3, "bond");
    ExcludedVolumeForceFieldGenerator::index_pairs_type contact_ignore_list =
        topology.ignore_list_within_edge(1, "contact");
    ignore_list.insert(ignore_list.end(),
          contact_ignore_list.begin(), contact_ignore_list.end());
    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    return ExcludedVolumeForceFieldGenerator(
            eps, 2.0/*cutoff ratio*/, radius_vec, ignore_list, use_periodic);
}

