#ifndef OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP

#include <regex>
#include <OpenMM.h>
#include "src/Simulator.hpp"
#include "src/SystemGenerator.hpp"
#include "src/Topology.hpp"
#include "ReadGenesisForceFieldGenerator.hpp"

std::map<std::string, std::map<std::string, std::string>> read_inp_file(const std::string& inp_file_name)
{
    std::cerr << "reading " << inp_file_name << "..." << std::endl;
    std::ifstream ifs(inp_file_name);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + inp_file_name);
    }

    std::map<std::string, std::map<std::string, std::string>> inp_data;
    std::string line, section_name;
    while(std::getline(ifs, line))
    {
        std::smatch result;
        if(section_name.empty())
        {
            if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
            {
                section_name = result.str(1);
                inp_data.insert(std::make_pair(section_name, std::map<std::string, std::string>()));
            }
        }
        else
        {
            std::map<std::string, std::string>& table_contents = inp_data.at(section_name);
            if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
            {
                section_name = result.str(1);
                inp_data.insert(std::make_pair(section_name, std::map<std::string, std::string>()));
            }
            else if(!std::regex_match(line, std::regex("^\\s*$"))) // non-comment line
            {
                std::size_t hash_position = line.find_first_of("#");
                if(hash_position == std::string::npos)
                {
                    hash_position = line.length();
                }
                const std::string noncomment_part = line.substr(0, hash_position);
                std::regex_match(noncomment_part, result,
                                 std::regex("\\s*(\\S+)\\s*=\\s*(\\S+)\\s*"));
                table_contents.insert(std::make_pair(result.str(1), result.str(2)));
            }
        }
    }

    ifs.close();

    return inp_data;
}

std::vector<std::string> preprocess_top_file(
        const std::string& top_file_name, const std::string& file_path)
{
    std::ifstream ifs(file_path + top_file_name);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + top_file_name);
    }

    std::vector<std::string> file_contents;
    std::smatch result;
    std::string line;
    while(std::getline(ifs, line))
    {
        // TODO: add case for define macro
        if(std::regex_match(line, std::regex("^\\s*$")) || std::regex_match(line, std::regex("^;.*")))
        {
            continue; // skip blank or comment line
        }
        else if(std::regex_match(line, result, std::regex("^#include \"(\\S*)\"\\s*")))
        {
            // expand include line
            const std::string include_file_name = result.str(1);
            std::cerr << "    reading include file " << file_path
                      << include_file_name << "..." << std::endl;
            std::vector<std::string> include_file = preprocess_top_file(include_file_name, file_path);
            for(std::string& line : include_file)
            {
                file_contents.push_back(line);
            }
        }
        else
        {
            file_contents.push_back(line);
        }
    }

    return file_contents;
}

std::map<std::string, std::vector<std::string>>
read_top_file(const std::string& topfile_name, const std::string& file_path)
{
    std::cerr << "reading " << file_path << topfile_name << "..." << std::endl;

    const std::vector<std::string>& file_contents = preprocess_top_file(topfile_name, file_path);

    std::map<std::string, std::vector<std::string>> top_data;
    std::string section_name, next_section_name;
    std::size_t line_idx(0), total_line_num(file_contents.size());
    while(line_idx < total_line_num)
    {
        std::smatch result;
        std::string line = file_contents[line_idx];
        if(!next_section_name.empty())
        {
            section_name = next_section_name;
            top_data.insert(std::make_pair(section_name, std::vector<std::string>()));
            next_section_name.clear();
        }
        else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
        {
            section_name = result.str(1);
            top_data.insert(std::make_pair(section_name, std::vector<std::string>()));
            line = file_contents[++line_idx];
        }
        else
        {
            throw std::runtime_error("[error] the line \"" + line +
                    "\" does not belong to any section in " + topfile_name + " file.");
        }

        // read contents of each section
        std::vector<std::string>& section_contents = top_data.at(section_name);
        while(line_idx < total_line_num)
        {
            if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
            {
                next_section_name = result.str(1);
                ++line_idx;
                break;
            }
            else
            {
                section_contents.push_back(line);
                line = file_contents[++line_idx];
            }
        }
        section_name = next_section_name;
    }

    return top_data;
}

std::vector<OpenMM::Vec3> read_genesis_initial_conf(
        const std::string& input_gro, const std::string& file_path)
{
    std::ifstream ifs(file_path + input_gro);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + input_gro);
    }

    std::string line;
    std::getline(ifs, line); // skip first line
    std::getline(ifs, line); // get second line
    const std::size_t system_size = std::stoi(line);
    std::vector<OpenMM::Vec3> initPosInNm(system_size);
    for(std::size_t idx=0; idx<system_size; ++idx)
    {
        std::getline(ifs, line);
        initPosInNm[idx] = OpenMM::Vec3(std::stof(line.substr(20, 9)),
                                        std::stof(line.substr(29, 9)),
                                        std::stof(line.substr(38, 9))); // location, nm
    }

    return initPosInNm;
}

Simulator make_simulator_from_genesis_inputs(
        const std::map<std::string, std::map<std::string, std::string>>& inpfile_data,
        const std::string& topfile_name, const std::string& grofile_name,
        const std::string& file_path)
{
    // read [INPUT] section
    const std::map<std::string, std::vector<std::string>>& top_data =
        read_top_file(topfile_name, file_path);

    if(top_data.find("atoms") == top_data.end())
    {
        throw std::runtime_error(
                "[error] There is no atoms section. Top file must contain one atoms section.");
    }
    std::vector<double>      mass_vec;
    std::vector<std::string> res_name_vec, atom_name_vec;
    for(auto& atoms_line : top_data.at("atoms"))
    {
        mass_vec.push_back(std::stof(atoms_line.substr(49, 9))); // amu
        res_name_vec .push_back(atoms_line.substr(27, 3));
        atom_name_vec.push_back(atoms_line.substr(31, 4));
    }
    SystemGenerator system_gen(mass_vec);

    // read [BOUNDARY] section
    const std::map<std::string, std::string> boundary_section = inpfile_data.at("BOUNDARY");
    const std::string type                   = boundary_section.at("type");
    const bool        use_periodic           = (type == "PBC");
    if(use_periodic)
    {
        std::cerr << "    boundary type is periodic boundary condition" << std::endl;
    }

    std::cerr << "reading forcefield tables..." << std::endl;
    Topology topology(top_data.at("atoms").size());
    if(top_data.find("bonds") != top_data.end())
    {
        HarmonicBondForceFieldGenerator ff_gen =
            read_genesis_harmonic_bond_ff_generator(top_data.at("bonds"), topology, use_periodic);
        if(ff_gen.indices().size() != 0)
        {
            system_gen.add_ff_generator(
                    std::make_unique<HarmonicBondForceFieldGenerator>(ff_gen));
        }
        else
        {
            std::cerr << "        -> skip this forcefield generation" << std::endl;
        }
    }

    if(top_data.find("angles") != top_data.end())
    {
        // make force field generator for AICG2+ angle
        GaussianBondForceFieldGenerator aicg_ff_gen =
            read_genesis_gaussian_bond_ff_generator(top_data.at("angles"), use_periodic);
        if(aicg_ff_gen.indices().size() != 0)
        {
            system_gen.add_ff_generator(
                    std::make_unique<GaussianBondForceFieldGenerator>(aicg_ff_gen));
        }
        else
        {
            std::cerr << "        -> skip this forcefield generation" << std::endl;
        }

        // make force field generator for Flexible Local angle
        for(const auto& aa_type_table : Constant::fla_spline_table)
        {
            FlexibleLocalAngleForceFieldGenerator fla_ff_gen =
                read_genesis_flexible_local_angle_ff_generator(top_data.at("angles"),
                                                               top_data.at("atoms"),
                                                               aa_type_table.first,
                                                               use_periodic);
            if(fla_ff_gen.indices().size() != 0)
            {
                system_gen.add_ff_generator(
                        std::make_unique<FlexibleLocalAngleForceFieldGenerator>(
                            fla_ff_gen));
            }
            else
            {
                std::cerr << "        -> skip this forcefield generation" << std::endl;
            }
        }
    }

    if(top_data.find("dihedrals") != top_data.end())
    {
        // make force field generator for AICG2+ dihedral
        GaussianDihedralForceFieldGenerator aicg_ff_gen =
            read_genesis_gaussian_dihedral_ff_generator(top_data.at("dihedrals"), use_periodic);
        if(aicg_ff_gen.indices().size() != 0)
        {
            system_gen.add_ff_generator(
                    std::make_unique<GaussianDihedralForceFieldGenerator>(aicg_ff_gen));
        }
        else
        {
            std::cerr << "        -> skip this forcefield generation" << std::endl;
        }

        // make force field generator for Flexible Local dihedral
        for(const auto& aa_type_pair_table : Constant::fld_fourier_table)
        {
             FlexibleLocalDihedralForceFieldGenerator fld_ff_gen =
                 read_genesis_flexible_local_dihedral_ff_generator(top_data.at("dihedrals"),
                                                                   top_data.at("atoms"),
                                                                   aa_type_pair_table.first,
                                                                   use_periodic);
             if(fld_ff_gen.indices().size() != 0)
             {
                 system_gen.add_ff_generator(
                     std::make_unique<FlexibleLocalDihedralForceFieldGenerator>(
                         fld_ff_gen));
             }
             else
             {
                 std::cerr << "        -> skip this forcefield generation." << std::endl;
             }
        }
    }

    if(top_data.find("pairs") != top_data.end())
    {
        GoContactForceFieldGenerator ff_gen =
            read_genesis_go_contact_ff_generator(top_data.at("pairs"), topology, use_periodic);
        if(ff_gen.indices().size() != 0)
        {
            system_gen.add_ff_generator(
                    std::make_unique<GoContactForceFieldGenerator>(ff_gen));
        }
        else
        {
            std::cerr << "        -> skip this forcefield generation" << std::endl;
        }
    }

    if(top_data.find("atomtypes") != top_data.end())
    {
        if(top_data.find("moleculetype") == top_data.end())
        {
            throw std::runtime_error("[error] There is no `[ moleculetype ]` section. Genesis input mode needs this section.");
        }
        const std::vector<std::string>& moleculetype_data = top_data.at("moleculetype");
        if(moleculetype_data.size() != 1)
        {
            throw std::runtime_error("[error] `[moleculetype]` section have " +
                   std::to_string(moleculetype_data.size()) +
                   " lines. Genesis input mode intends this section have only 1 line." );
        }

        const std::size_t ignore_particle_within_bond =
            std::stoi(moleculetype_data[0].substr(17, 6));
        ExcludedVolumeForceFieldGenerator ff_gen =
            read_genesis_exv_ff_generator(top_data.at("atomtypes"),
                top_data.at("atoms"), topology, use_periodic, ignore_particle_within_bond);
        system_gen.add_ff_generator(
                std::make_unique<ExcludedVolumeForceFieldGenerator>(ff_gen));
    }
    else
    {
        throw std::runtime_error("[error] There is no `[ atomtypes ]` section. Genesis input mode needs this section.");
    }

    const std::vector<OpenMM::Vec3> initial_position_in_nm(
            read_genesis_initial_conf(grofile_name, file_path));

    // read [OUTPUT] section
    const std::map<std::string, std::string>& output_section = inpfile_data.at("OUTPUT");
    const std::string& dcdfile_name = output_section.at("dcdfile");
    const std::string& pdbfile_name = output_section.at("pdbfile");
    std::size_t file_suffix_from = dcdfile_name.rfind(".");
    if(file_suffix_from == std::string::npos)
    {
        file_suffix_from = dcdfile_name.size();
    }
    const std::string file_prefix  = dcdfile_name.substr(0, file_suffix_from);
    const std::string enefile_name = file_prefix + ".ene";

    // read [DYNAMICS] section
    const std::map<std::string, std::string> dynamics_section = inpfile_data.at("DYNAMICS");
    const std::size_t nsteps        = std::stoi(dynamics_section.at("nsteps"));
    const double      timestep      = std::stof(dynamics_section.at("timestep")); // ps
    const std::size_t crdout_period = std::stoi(dynamics_section.at("crdout_period"));
    std::size_t seed = 0;
    if(dynamics_section.find("iseed") != dynamics_section.end())
    {
        seed = std::stoi(dynamics_section.at("iseed"));
    }

    // read [EMSEMBLE] section
    const std::map<std::string, std::string> emsemble_section = inpfile_data.at("ENSEMBLE");
    const double temperature = std::stof(emsemble_section.at("temperature"));
    const double gamma_t     = std::stof(emsemble_section.at("gamma_t")); // GENESIS gamma_t unit is ps^-1

    // setup OpenMM integrator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for
    // different molecules. However, cafemol use different gamma depends on the
    // mass of each particle, and the product of mass and gamma is constant in
    // there. We need to implement new LangevinIntegrator which can use different
    // gamma for different particles.
    auto integrator = OpenMM::LangevinIntegrator(temperature,
                                                 gamma_t/*friction coef ps^-1*/,
                                                 timestep/* delta t */);
    integrator.setRandomNumberSeed(seed);

    // construct observers
    std::vector<std::unique_ptr<ObserverBase>> observers;
    observers.push_back(std::make_unique<DCDObserver>(
                file_prefix, nsteps, crdout_period, timestep, use_periodic));
    observers.push_back(std::make_unique<EnergyObserver>(file_prefix, system_gen));

    // dump initial configuration to pdb file
    std::cerr << "    output initial state file : " << pdbfile_name << std::endl;
    Utility::clear_file(pdbfile_name);
    std::ofstream ofs(pdbfile_name, std::ios::app);
    ofs << "MODEL     0" << std::endl;
    for(std::size_t idx=0; idx<initial_position_in_nm.size(); ++idx)
    {
        ofs << std::setprecision(3);
        ofs << "ATOM  "                     //                          1-6
            << std::setw(5) << idx+1 << " " // atom serial number       7-12
            << atom_name_vec[idx] << " "    // atom name               13-17
            << res_name_vec [idx] << "  "   // residue name            18-22
            << "   1    ";                  // residue sequence number 23-30
        ofs << std::setw(8) << std::fixed
            << std::setw(8) << std::fixed << initial_position_in_nm[idx][0]*OpenMM::AngstromsPerNm
            << std::setw(8) << std::fixed << initial_position_in_nm[idx][1]*OpenMM::AngstromsPerNm
            << std::setw(8) << std::fixed << initial_position_in_nm[idx][2]*OpenMM::AngstromsPerNm
            << "  1.00  0.00" << std::endl;
    }
    ofs << "ENDMDL" << std::endl; // end of frame
    ofs.close();

    return Simulator(system_gen, integrator,
                     initial_position_in_nm, nsteps, crdout_period, observers);
}

Simulator read_genesis_input(const std::string& inp_file_name)
{
    std::size_t file_path_len = inp_file_name.rfind("/")+1;
    if(file_path_len == std::string::npos)
    {
        file_path_len = 0;
    }
    const std::string file_path = inp_file_name.substr(0, file_path_len);

    const std::map<std::string, std::map<std::string, std::string>>& inp_data =
        read_inp_file(inp_file_name);

    const std::map<std::string, std::string>& input_section = inp_data.at("INPUT");
    const std::string& topfile_name = input_section.at("grotopfile");
    const std::string& grofile_name = input_section.at("grocrdfile");

    return make_simulator_from_genesis_inputs(inp_data, topfile_name, grofile_name, file_path);
}

#endif // OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
