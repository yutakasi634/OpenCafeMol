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

    struct bond_param
    {
        const std::pair<std::size_t, std::size_t> indices_;
        const double                              equilibrium_value_;
        const double                              force_coef_;

        bond_param(const std::pair<std::size_t, std::size_t> idxs,
                const double eq_val, const double f_coef)
            : indices_(idxs), equilibrium_value_(eq_val), force_coef_(f_coef)
        {}
    };

    std::vector<bond_param> type1_bond_params; // harmonic_bond
    std::vector<bond_param> type21_bond_params;

    for(const auto& bond_line : bonds_data)
    {
        const std::size_t idx_i  = std::stoi(bond_line.substr( 0, 10)) - 1;
        const std::size_t idx_j  = std::stoi(bond_line.substr(10, 10)) - 1;
        const std::size_t f_type = std::stoi(bond_line.substr(20,  5));
        const double      eq     = std::stod(bond_line.substr(25, 42));
        const double      coef   = std::stod(bond_line.substr(43, 60)); // KJ/(mol nm^2)
        if(f_type == 1)
        {
            type1_bond_params.push_back(
                    bond_param(std::make_pair(idx_i, idx_j), eq, coef * 0.5));
        }
        else if(f_type == 21)
        {
            type21_bond_params.push_back(
                    bond_param(std::make_pair(idx_i, idx_j), eq, coef));
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
            indices_vec.emplace_back(param.indices_);
            eqs        .emplace_back(param.equilibrium_value_);
            coefs      .emplace_back(param.force_coef_);
        }
        topology.add_edges(indices_vec, "bond");

        std::cerr << "    BondLength    : Harmonic (" <<
                  indices_vec.size() << " found)" << std::endl;

        ff_gen_ptrs.push_back(
                std::make_unique<HarmonicBondForceFieldGenerator>(
                    indices_vec, eqs, coefs, use_periodic));
    }

    if(type21_bond_params.size() != 0)
    {
        std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
        std::vector<double>                              eqs;
        std::vector<double>                              coefs;
        for(const auto& param: type21_bond_params)
        {
            indices_vec.emplace_back(param.indices_);
            eqs        .emplace_back(param.equilibrium_value_);
            coefs      .emplace_back(param.force_coef_);
        }
        topology.add_edges(indices_vec, "bond");

        std::cerr << "    BondLength    : 3SPN2 (" <<
                  indices_vec.size() << " found)" << std::endl;

        ff_gen_ptrs.push_back(
                std::make_unique<ThreeSPN2BondForceFieldGenerator>(
                    indices_vec, coefs, coefs, eqs, use_periodic));
    }

    return ff_gen_ptrs;
}

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_angles_section(
        const std::vector<std::string>& angles_data,
        const bool use_periodic, const std::vector<std::string>& res_name_vec)
{
    struct angle_param
    {
        const std::array<std::size_t, 3> indices_;
        const std::string                param_str_;

        angle_param(const std::array<std::size_t, 3> indices,
                const std::string param_str)
            : indices_(indices), param_str_(param_str)
        {}
    };

    std::vector<angle_param> type1_angle_params;  // harmonic_angle
    // type21 angle potential is represented by gaussian bond potential internally.
    // this is not a mistake.
    std::vector<angle_param> type21_angle_params; // gaussian_bond
    std::vector<angle_param> type22_angle_params; // flp_angle

    for(const auto& angle_line : angles_data)
    {
        const std::array<std::size_t, 3> indices = {
            std::stoul(angle_line.substr( 0, 10)) - 1,
            std::stoul(angle_line.substr(10, 10)) - 1,
            std::stoul(angle_line.substr(20, 10)) - 1
        };
        const std::size_t f_type    = std::stoul(angle_line.substr(30,  5));
        const std::string param_str = angle_line.substr(35);

        if(f_type == 1)
        {
            type1_angle_params.push_back(angle_param(indices, param_str));
        }
        else if(f_type == 21)
        {
            type21_angle_params.push_back(angle_param(indices, param_str));
        }
        else if(f_type == 22)
        {
            type22_angle_params.push_back(angle_param(indices, param_str));
        }
        else
        {
            throw std::runtime_error(
                    "[error] unexpected angle type " + std::to_string(f_type) +
                    " was specified in `[ angles ]` section."
                    " expected value is 1, 21 or 22.");
        }
    }

    std::vector<std::unique_ptr<ForceFieldGeneratorBase>> ff_gen_ptrs{};

    if(type1_angle_params.size() != 0)
    {
        std::vector<std::array<std::size_t, 3>> indices_vec;
        std::vector<double>                     thetas; // rad
        std::vector<double>                     ks;     // kJ/(mol rad^2)

        for(const auto& angle_param : type1_angle_params)
        {
            indices_vec.push_back(angle_param.indices_);
            const double theta_degree =
                std::stod(angle_param.param_str_.substr( 0, 15));  // °
            thetas.push_back(theta_degree / 180.0 * Constant::pi); // rad
            ks    .push_back(std::stod(angle_param.param_str_.substr(15, 15)));
        }

        std::cerr << "    BondAngle     : Harmonic (" << indices_vec.size()
                  << " found)" << std::endl;

        ff_gen_ptrs.push_back(
                std::make_unique<HarmonicAngleForceFieldGenerator>(
                    indices_vec, thetas, ks, use_periodic));
    }

    if(type21_angle_params.size() != 0)
    {
        std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
        std::vector<double>                              r0s; // nm
        std::vector<double>                              eps; // kJ/mol
        std::vector<double>                              ws;  // nm

        for(const auto& angle_param : type21_angle_params)
        {
            indices_vec.push_back(std::make_pair(
                        angle_param.indices_[0], angle_param.indices_[2]));
            const std::string& param_str = angle_param.param_str_;
            r0s.push_back(std::stod(param_str.substr( 0, 15)));
            eps.push_back(std::stod(param_str.substr(15, 15)));
            ws .push_back(std::stod(param_str.substr(30, 15)));
        }

        std::cerr << "    BondLength    : Gaussian (" << indices_vec.size()
                  << " found)" << std::endl;

        ff_gen_ptrs.push_back(
                std::make_unique<GaussianBondForceFieldGenerator>(
                    indices_vec, r0s, eps, ws, use_periodic));
    }

    if(type22_angle_params.size() != 0)
    {
        for(const auto& aa_type_table : Constant::fla_spline_table)
        {
            const std::string& aa_type = aa_type_table.first;

            std::vector<std::array<std::size_t, 3>> indices_vec;
            for(const auto& angle_param : type22_angle_params)
            {
                const std::size_t idx_j = angle_param.indices_[1];
                if(res_name_vec[idx_j] == aa_type)
                {
                    indices_vec.push_back(angle_param.indices_);
                }
            }

            std::cerr << "    BondAngle     : FlexibleLocalAngle - "
                      << aa_type << " (" << indices_vec.size() << " found)"
                      << std::endl;

            if(indices_vec.size() == 0)
            {
                std::cerr << "        -> skip this forcefield generation" << std::endl;
                continue;
            }

            // the unit of spline table is kcal/mol and OpenMM eneryg unit is kJ/mol.
            // this unit conversion is performed through ks;
            std::vector<double> ks =
                std::vector<double>(indices_vec.size(), OpenMM::KJPerKcal);
            ff_gen_ptrs.push_back(
                    std::make_unique<FlexibleLocalAngleForceFieldGenerator>(
                        indices_vec, ks,
                        Constant::fla_spline_table.at(aa_type), aa_type, use_periodic));
        }
    }

    return ff_gen_ptrs;
}

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_dihedrals_section(
        const std::vector<std::string>& dihedrals_data,
        const bool use_periodic, const std::vector<std::string>& res_name_vec)
{
    struct dihedral_param
    {
        const std::array<std::size_t, 4> indices_;
        const std::string                param_str_;

        dihedral_param(const std::array<std::size_t, 4> indices,
                const std::string param_str)
            : indices_(indices), param_str_(param_str)
        {}
    };

    // genesis use safe type potentials instead of original potentials.
    // this software temporarily processes these safe potentials as the original one.
    // TODO: we need to adjust safe potential case.
    std::vector<dihedral_param> type1_dihedral_params;  // cosine dihedral
    std::vector<dihedral_param> type21_dihedral_params; // gaussian dihedral
    std::vector<dihedral_param> type22_dihedral_params; // flp dihedral

    for(const auto& dihedral_line : dihedrals_data)
    {
        const std::array<std::size_t, 4> indices = {
            std::stoul(dihedral_line.substr( 0, 10)) - 1,
            std::stoul(dihedral_line.substr(10, 10)) - 1,
            std::stoul(dihedral_line.substr(20, 10)) - 1,
            std::stoul(dihedral_line.substr(30, 10)) - 1
        };
        const std::size_t f_type    = std::stoul(dihedral_line.substr(40, 5));
        const std::string param_str = dihedral_line.substr(45);

        if(f_type == 1 || f_type == 32)
        {
            type1_dihedral_params.push_back(dihedral_param(indices, param_str));
        }
        else if(f_type == 21 || f_type == 41)
        {
            type21_dihedral_params.push_back(dihedral_param(indices, param_str));
        }
        else if(f_type == 22 || f_type == 52)
        {
            type22_dihedral_params.push_back(dihedral_param(indices, param_str));
        }
        else
        {
            throw std::runtime_error(
                    "[error] unexpected angle type " + std::to_string(f_type) +
                    " was specified in `[ dihedrals ]` section."
                    " expected value is 1(32), 21(41) or 22(52).");
        }
    }

    std::vector<std::unique_ptr<ForceFieldGeneratorBase>> ff_gen_ptrs{};

    if(type1_dihedral_params.size() != 0)
    {
        std::vector<std::array<std::size_t, 4>> indices_vec;
        std::vector<double>                     phi0s; // rad
        std::vector<double>                     ks;    // kJ/mol
        std::vector<double>                     ns;    // dimensionless

        for(const auto& dihedral_param : type1_dihedral_params)
        {
            indices_vec.push_back(dihedral_param.indices_);
            const double phi0_degree =
                std::stod(dihedral_param.param_str_.substr(0, 15)); // °
            phi0s.push_back(phi0_degree / 180.0 * Constant::pi);    // rad
            ks   .push_back(std::stod(dihedral_param.param_str_.substr(15, 15)));
            ns   .push_back(std::stoi(dihedral_param.param_str_.substr(30, 15)));
        }

        std::cerr << "    DihedralAngle : Cosine (" << indices_vec.size()
                  << " found)" << std::endl;

        ff_gen_ptrs.push_back(std::make_unique<CosineDihedralForceFieldGenerator>(
                    indices_vec, ks, phi0s, ns, use_periodic));
    }

    if(type21_dihedral_params.size() != 0)
    {
        std::vector<std::array<std::size_t, 4>> indices_vec;
        std::vector<double>                     phi0s;  // rad
        std::vector<double>                     eps;    // kJ/mol
        std::vector<double>                     sigmas; // rad

        for(const auto& dihedral_param : type21_dihedral_params)
        {
            indices_vec.push_back(dihedral_param.indices_);
            const double phi0_degree =
                std::stod(dihedral_param.param_str_.substr(0, 15)); // °
            phi0s .push_back(phi0_degree / 180.0 * Constant::pi);    // rad
            eps   .push_back(std::stod(dihedral_param.param_str_.substr(15, 15)));
            sigmas.push_back(std::stod(dihedral_param.param_str_.substr(30, 15)));
        }

        std::cerr << "    DihedralAngle : Gaussian (" << indices_vec.size() 
                  << " found)" << std::endl;


        ff_gen_ptrs.push_back(std::make_unique<GaussianDihedralForceFieldGenerator>(
                    indices_vec, eps, phi0s, sigmas, use_periodic));
    }

    if(type22_dihedral_params.size() != 0)
    {
        for(const auto& aa_type_pair_table : Constant::fld_fourier_table)
        {
            const std::pair<std::string, std::string> aa_type_pair =
                aa_type_pair_table.first;

            std::vector<std::array<std::size_t, 4>> indices_vec;
            for(const auto& dihedral_param : type22_dihedral_params)
            {
                const std::array<std::size_t, 4>& indices = dihedral_param.indices_;
                if(aa_type_pair.first == "R1")
                {
                    if(res_name_vec[indices[2]] == "GLY")
                    {
                        indices_vec.push_back(indices);
                    }
                }
                else if(aa_type_pair.first == "R2")
                {
                    if(res_name_vec[indices[1]] != "GLY" &&
                        res_name_vec[indices[2]] == "PRO")
                    {
                        indices_vec.push_back(indices);
                    }
                }
                else if(aa_type_pair.first == "GLY" && aa_type_pair.second == "PRO")
                {
                    if(res_name_vec[indices[1]] == "GLY" &&
                            res_name_vec[indices[2]] == "PRO")
                    {
                        indices_vec.push_back(indices);
                    }
                }
                else
                {
                    if(res_name_vec[indices[2]] == aa_type_pair.first)
                    {
                        indices_vec.push_back(indices);
                    }
                }
            }

            std::cerr << "    DihedralAngle : FlexibleLocalDihedral - "
                      << aa_type_pair.first << "-" << aa_type_pair.second
                      << " (" << indices_vec.size() << " found)" << std::endl;

            if(indices_vec.size() == 0)
            {
                std::cerr << "        -> skip this forcefield generation" << std::endl;
                continue;
            }

            // the unit of fourier table is kcal/mol and OpenMM eneryg unit is kJ/mol.
            // this unit conversion is performed through ks;
            std::vector<double> ks =
                std::vector<double>(indices_vec.size(), OpenMM::KJPerKcal);
            ff_gen_ptrs.push_back(
                    std::make_unique<FlexibleLocalDihedralForceFieldGenerator>(
                        indices_vec, ks, Constant::fld_fourier_table.at(aa_type_pair),
                        aa_type_pair.first + "-" + aa_type_pair.second, use_periodic));
        }
    }

    return ff_gen_ptrs;
}

GoContactForceFieldGenerator
read_genesis_pairs_section(
        const std::vector<std::string>& pairs_data, Topology& topology,
        const bool use_periodic)
{
    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              ks;
    std::vector<double>                              r0s;

    for(auto& pairs_line : pairs_data)
    {
        const std::size_t idx_i  = std::stoul(pairs_line.substr(0,  10)) - 1;
        const std::size_t idx_k  = std::stoul(pairs_line.substr(10, 10)) - 1;
        const std::size_t f_type = std::stoul(pairs_line.substr(20, 10));
        const double      r0     = std::stod(pairs_line.substr(30, 15)); // nm
        const double      k      = std::stod(pairs_line.substr(45, 15)); // KJ/mol

        if(f_type != 2)
        {
            throw std::runtime_error(
                    "[error] unexpected pair type " + std::to_string(f_type) +
                    " was specified in `[ pairs ]` section."
                    " expected value is only 2.");
        }

        indices_vec.push_back(std::make_pair(idx_i, idx_k));
        ks         .push_back(k);
        r0s        .push_back(r0);
    }
    topology.add_edges(indices_vec, "contact");

    std::cerr << "    BondLength    : GoContact (" << indices_vec.size() << " found)" << std::endl;
    return GoContactForceFieldGenerator(indices_vec, ks, r0s, use_periodic);
}

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_basestacking_ff_generator(
        const std::vector<std::string>& basestacktypes_data,
        const std::vector<std::string>& atoms_data, Topology& topology,
        const bool use_periodic)
{
    // read `[ basestacktypes ]` section
}

std::vector<std::unique_ptr<ForceFieldGeneratorBase>>
read_genesis_exv_ff_generators(
        const std::map<std::string, std::vector<std::string>>& top_data,
        const std::vector<std::string>& particle_type_vec,
        Topology& topology, const bool use_periodic,
        const std::size_t ignore_particle_within_bond)
{
    std::vector<std::unique_ptr<ForceFieldGeneratorBase>> ff_gen_ptrs{};

    // read `[ atomtypes ]` section
    const std::vector<std::string>& atomtypes_data = top_data.at("atomtypes");
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
    for(auto& particle_type : particle_type_vec)
    {
        radius_vec.push_back(name_rmin_map.at(particle_type));
    }

    std::vector<std::optional<std::string>> group_vec(topology.size(), std::nullopt);
    std::vector<std::pair<std::string, std::string>> exv_ignore_group_pairs;
    if(top_data.find("cgdnaexvtypes") != top_data.end())
    {
        const std::vector<std::string>& cgdnas_data = top_data.at("cgdnaexvtypes");
        std::map<std::string, double> name_sigma_map;
        for(auto& cgdna_line : cgdnas_data)
        {
            const std::size_t f_type = std::stoul(cgdna_line.substr(6, 5));
            if(f_type != 1)
            {
                throw std::runtime_error(
                        "[error] unexpected func type " + std::to_string(f_type) +
                        " was specified in `[ cgdnaexvtypes ]` section. "
                        "expected value is 1.");
            }
            // make dna atom and func, sigma list
            const std::string name  = Utility::erase_space(cgdna_line.substr(0, 6));
            const double      sigma = std::stod(cgdna_line.substr(12, 10));
            name_sigma_map.insert(std::make_pair(name, sigma));
        }

        // setup ignore group for EXVFFGen
        for(std::size_t p_idx=0; p_idx < particle_type_vec.size(); ++p_idx)
        {
            const std::string& particle_type(particle_type_vec[p_idx]);
            if(name_sigma_map.find(particle_type) != name_sigma_map.end())
            {
                group_vec[p_idx] = ("DNA");
            }
        }
        exv_ignore_group_pairs.push_back(std::make_pair("DNA", "DNA"));

        // setup 3SPN2EXFVFFGen
        std::size_t exv_3spn2_particle_count(0);
        std::vector<std::optional<double>> exv_3spn2_radius_vec;
        for(const auto& particle_type : particle_type_vec)
        {
            if(name_sigma_map.find(particle_type) != name_sigma_map.end())
            {
                exv_3spn2_radius_vec.push_back(name_sigma_map.at(particle_type));
                ++exv_3spn2_particle_count;
            }
            else
            {
                exv_3spn2_radius_vec.push_back(std::nullopt);
            }
        }

        if(exv_3spn2_particle_count != 0)
        {
            std::vector<std::pair<std::size_t, std::size_t>> ignore_list =
                topology.ignore_list_within_edge(1, "bond");
            std::vector<std::pair<std::size_t, std::size_t>> nuc_ignore_list =
                topology.ignore_list_within_edge(1, "nucleotide");
            ignore_list.insert(ignore_list.end(),
                    nuc_ignore_list.begin(), nuc_ignore_list.end());

            // for exclude base pairing interaction pair
            for(std::size_t idx=0; idx<particle_type_vec.size(); ++idx)
            {
                for(std::size_t jdx=idx+1; jdx<particle_type_vec.size(); ++jdx)
                {
                    const auto& type_i = particle_type_vec[idx];
                    const auto& type_j = particle_type_vec[jdx];

                    // complement
                    if ((type_i == "DA" && type_j == "DT") ||
                        (type_i == "DT" && type_j == "DA") ||
                        (type_i == "DG" && type_j == "DC") ||
                        (type_i == "DC" && type_j == "DG"))
                    {
                        ignore_list.push_back(std::make_pair(idx, jdx));
                    }
                }
            }

            std::sort(ignore_list.begin(), ignore_list.end());
            const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
            ignore_list.erase(result, ignore_list.end());

            std::cerr << "    Global        : ExcludedVolume 3SPN2 ("
                      << exv_3spn2_particle_count << " found)" << std::endl;

            const auto eps = ThreeSPN2ExcludedVolumePotentialParameter::epsilon;
            ff_gen_ptrs.push_back(
                    std::make_unique<ThreeSPN2ExcludedVolumeForceFieldGenerator>(
                        eps, 2.0/*cutoff ratio*/, radius_vec,
                        ignore_list, use_periodic));
        }
    }

    std::cerr << "    Global        : ExcludedVolume (" << radius_vec.size()
              << " found)"                              << std::endl;
    // generate ignore list
    ExcludedVolumeForceFieldGenerator::index_pairs_type ignore_list =
        topology.ignore_list_within_edge(ignore_particle_within_bond, "bond");
    ExcludedVolumeForceFieldGenerator::index_pairs_type contact_ignore_list =
        topology.ignore_list_within_edge(1, "contact");
    ignore_list.insert(ignore_list.end(),
          contact_ignore_list.begin(), contact_ignore_list.end());
    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    ff_gen_ptrs.push_back(std::make_unique<ExcludedVolumeForceFieldGenerator>(
            eps, 2.0/*cutoff ratio*/, radius_vec, ignore_list, use_periodic,
            exv_ignore_group_pairs, group_vec));



    return ff_gen_ptrs;
}

