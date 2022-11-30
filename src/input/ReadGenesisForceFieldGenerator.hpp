#ifndef OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "src/forcefield/GaussianBondForceFieldGenerator.hpp"
#include "src/forcefield/GoContactForceFieldGenerator.hpp"
#include "src/forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "src/forcefield/GaussianDihedralForceFieldGenerator.hpp"

const HarmonicBondForceFieldGenerator
read_genesis_harmonic_bond_ff_generator(const std::vector<std::string>& bonds_data)
{
        std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
        std::vector<double>                              v0s;
        std::vector<double>                              ks;

        for(auto& bonds_line : bonds_data)
        {
            const std::size_t idx_i = std::stoi(bonds_line.substr( 0, 10)) - 1;
            const std::size_t idx_j = std::stoi(bonds_line.substr(10, 10)) - 1;
            const double      v0    = std::stof(bonds_line.substr(25, 42));
            const double      k     = std::stof(bonds_line.substr(43, 60)) * 0.5; // KJ/(mol nm^2)
            indices_vec.push_back(std::make_pair(idx_i, idx_j));
            v0s        .push_back(v0);
            ks         .push_back(k);
        }

        std::cerr << "    BondLength    : Harmonic (" <<
                      indices_vec.size() << " found)" << std::endl;
        return HarmonicBondForceFieldGenerator(indices_vec, v0s, ks);
}

const GaussianBondForceFieldGenerator
read_genesis_gaussian_bond_ff_generator(const std::vector<std::string>& angles_data)
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
            const double      v0    = std::stof(angles_line.substr(35, 15)); // nm
            const double      k     = std::stof(angles_line.substr(50, 15)); // KJ/mol
            const double      sigma = std::stof(angles_line.substr(65, 15)); // nm

            indices_vec.push_back(std::make_pair(idx_i, idx_k));
            v0s        .push_back(v0);
            ks         .push_back(k);
            sigmas     .push_back(sigma);
        }
    }

    std::cerr << "    BondLength    : Gaussian (" << indices_vec.size() << " found)" << std::endl;
    return GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas);
}

const GoContactForceFieldGenerator
read_genesis_go_contact_ff_generator(const std::vector<std::string>& pairs_data)
{
    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(auto& pairs_line : pairs_data)
    {
        const std::size_t idx_i = std::stoi(pairs_line.substr(0,  10)) - 1;
        const std::size_t idx_k = std::stoi(pairs_line.substr(10, 10)) - 1;
        const double      r0    = std::stof(pairs_line.substr(30, 15)); // nm
        const double      k     = std::stof(pairs_line.substr(45, 15)); // KJ/mol

        indices_vec.push_back(std::make_pair(idx_i, idx_k));
        ks         .push_back(k);
        r0s        .push_back(r0);
    }

    std::cerr << "    BondLength    : GoContact (" << indices_vec.size() << " found)" << std::endl;
    return GoContactForceFieldGenerator(indices_vec, ks, r0s);
}

const FlexibleLocalAngleForceFieldGenerator
read_genesis_flexible_local_angle_ff_generator(const std::vector<std::string>& angles_data,
        const std::vector<std::string>& atoms_data, const std::string& aa_type)
{
    std::vector<std::vector<std::size_t>> indices_vec;

    for(const auto& angles_line : angles_data)
    {
        if(std::stoi(angles_line.substr(30, 5)) == 22) // 22 means FLP-type angle function
        {
            const std::size_t idx_j = std::stoi(angles_line.substr(10, 10)) - 1;
            const std::string atom_type = erase_space(atoms_data[idx_j].substr(25, 5));
            if(atom_type == aa_type)
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
               indices_vec, std::vector<double>(indices_vec.size(), 1.0),
               Constant::fla_spline_table.at(aa_type), aa_type);
}


#endif // OPEN_AICG2_PLUS_READ_GENESIS_FORCE_FIELD_GENERATOR_HPP
