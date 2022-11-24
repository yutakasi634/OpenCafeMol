#ifndef OPEN_AICG2_PLUS_READ_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_INPUT_HPP

#include <memory>
#include <regex>
#include <OpenMM.h>
#include "src/Simulator.hpp"
#include "src/Topology.hpp"
#include "ReadForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::System> read_system(const toml::value& data)
{
    auto system_ptr = std::make_unique<OpenMM::System>();

    // read systems tables
    const auto& systems   = toml::find(data, "systems");
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::cerr << "generating system with size " << system_size << "..." << std::endl;
    for(std::size_t i=0; i<system_size; ++i)
    {
        // set mass
        const auto& p = particles.at(i);
        system_ptr->addParticle(toml::get<double>(find_either(p, "m", "mass"))); // amu
    }

    // for exclusion list of Excluded Volume
    std::vector<std::pair<std::size_t, std::size_t>> bonded_pairs;
    std::vector<std::pair<std::size_t, std::size_t>> contacted_pairs;

    std::cerr << "generating forcefields..." << std::endl;
    // read forcefields info
    const auto ff = toml::find(data, "forcefields").at(0);
    Topology topology(system_size);
    if(ff.contains("local"))
    {
        const auto& locals = toml::find(ff, "local").as_array();

        for(const auto& local_ff : locals)
        {
            const std::string interaction = toml::find<std::string>(local_ff, "interaction");
            const std::string potential   = toml::find<std::string>(local_ff, "potential");
            if(interaction == "BondLength" && potential == "Harmonic")
            {
                const auto ff_gen =
                    read_harmonic_bond_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondLength" && potential == "Gaussian")
            {
                const auto ff_gen =
                    read_gaussian_bond_ff_generator(local_ff);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondLength" && potential == "GoContact")
            {
                const auto ff_gen =
                    read_go_contact_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                for(const auto& [aa_type, spline_table] : Constant::fla_spline_table)
                {
                    const auto ff_gen =
                        read_flexible_local_angle_ff_generator(
                                local_ff, aa_type, spline_table);
                    system_ptr->addForce(ff_gen.generate().release());
                }
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                const auto ff_gen =
                    read_gaussian_dihedral_ff_generator(local_ff);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "DihedralAngle" && potential == "FlexibleLocalDihedral")
            {
                for(const auto& [aa_pair_type, fourier_table] : Constant::fld_fourier_table)
                {
                    const auto ff_gen =
                        read_flexible_local_dihedral_ff_generator(
                            local_ff, aa_pair_type, fourier_table);
                    system_ptr->addForce(ff_gen.generate().release());
                }
            }
        }
    }
    topology.make_molecule("bond");

    if(ff.contains("global"))
    {
        const auto& globals = toml::find(ff, "global").as_array();

        for(const auto& global_ff : globals)
        {
            const std::string potential = toml::find<std::string>(global_ff, "potential");
            if(potential == "ExcludedVolume")
            {
                const auto ff_gen =
                    read_excluded_volume_ff_generator(global_ff, system_size, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            if(potential == "DebyeHuckel")
            {
                const auto attr = toml::find(systems[0], "attributes");

                const auto temperature = toml::expect<double>(attr, "temperature");
                if(!temperature.is_ok())
                {
                    std::cerr << temperature.unwrap_err() << std::endl;
                }

                const auto ionic_strength = toml::expect<double>(attr, "ionic_strength");
                if(!ionic_strength.is_ok())
                {
                    std::cerr << ionic_strength.unwrap_err() << std::endl;
                }

                const auto ff_gen =
                    read_debye_huckel_ff_generator(global_ff, system_size,
                        ionic_strength.unwrap(), temperature.unwrap(), topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
        }
    }

    return system_ptr;
}

std::vector<OpenMM::Vec3> read_initial_conf(const toml::value& data)
{
    const auto& systems     = toml::find(data, "systems");
    const auto& particles   = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initPosInNm(system_size);
    for (std::size_t i=0; i<system_size; ++i)
    {
        // set position
        const auto& p = particles.at(i);
        std::array<double, 3> vec = {0.0, 0.0, 0.0};
        vec = toml::get<std::array<double, 3>>(find_either(p, "pos", "position"));
        initPosInNm[i] = OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                      vec[1]*OpenMM::NmPerAngstrom,
                                      vec[2]*OpenMM::NmPerAngstrom); // location, nm
    }

    return initPosInNm;
}

const std::string get_file_suffix(const std::string& filename)
{
    const std::size_t file_suffix_from = filename.rfind(".");
    if(file_suffix_from == std::string::npos)
    {
        throw std::runtime_error(
                "[error] There is no file extension in " + filename + "."
                " The file type can not be specified.");
    }
    const std::size_t file_suffix_len  = filename.length() - file_suffix_from;
    return filename.substr(file_suffix_from, file_suffix_len);
}

Simulator read_toml_input(const std::string& input_file_name)
{
    std::size_t file_path_len   = input_file_name.rfind("/")+1;
    if(file_path_len == std::string::npos)
    {
        file_path_len = 0;
    }

    const std::string file_suffix = get_file_suffix(input_file_name);
    if(file_suffix != ".toml")
    {
            throw std::runtime_error(
                    "[error] File suffix is `" + file_suffix + "`."
                    " Toml mode needs `.toml` file for input.");
    }
    const std::size_t file_prefix_from = input_file_name.rfind(".");
    const std::size_t file_prefix_len = file_prefix_from - file_path_len;
    const std::string file_path       = input_file_name.substr(0, file_path_len);
    const std::string file_prefix     = input_file_name.substr(file_path_len, file_prefix_len);

    // read input toml file
    std::cerr << "parsing " << input_file_name << "..." << std::endl;
    auto data = toml::parse(input_file_name);

    // read system table
    const auto&    systems     = toml::find(data, "systems");
    const auto&    attr        = toml::find(systems[0], "attributes");
    const auto&    temperature = toml::find<double>(attr, "temperature");
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_initial_conf(data));

    // read simulator table
    const auto&       simulator_table = toml::find(data, "simulator");
    const std::size_t total_step      = toml::find<std::size_t>(simulator_table, "total_step");
    const std::size_t save_step       = toml::find<std::size_t>(simulator_table, "save_step");
    const double      delta_t         = toml::find<double>(simulator_table, "delta_t");

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.2 ps^-1 temporary, correspond to
    // approximatry 0.01 in cafemol friction coefficient. We need to implement new
    // LangevinIntegrator which can use different gamma for different particles.
    return Simulator(std::move(read_system(data)),
               OpenMM::LangevinIntegrator(temperature,
                                          0.3/*friction coef ps^-1*/,
                                          delta_t*Constant::cafetime),
               initial_position_in_nm, total_step, save_step, Observer(file_prefix));
}

void read_genesis_input(
        const std::string& input_gro, const std::string& input_itp, const std::string& input_top)
{
    std::cerr << "gro file is " << input_gro << std::endl;
    std::cerr << "itp file is " << input_itp << std::endl;
    std::cerr << "top file is " << input_top << std::endl;

    // read input itp file
    std::cerr << "parsing " << input_itp << "..." << std::endl;
    std::ifstream ifs(input_itp);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + input_itp);
    }

    std::string line, next_section_name;
    while(std::getline(ifs, line))
    {
        std::smatch result;
        std::string section_name;
        if(!next_section_name.empty())
        {
            section_name = next_section_name;
            next_section_name.clear();
        }
        else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
        {
            section_name = result.str(1);
        }

        if(!section_name.empty())
        {
            std::cerr << section_name << " section detected" << std::endl;
            if(section_name == "moleculetype")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in moleculetype table " << line << std::endl;
                    }
                }
            }
            else if(section_name == "atoms")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in atom table " << line << std::endl;
                    }
                }
            }
            else if(section_name == "bonds")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in bonds table " << line << std::endl;
                    }
                }
            }
            else if(section_name == "angles")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in angles table " << line << std::endl;
                    }
                }
            }
            else if(section_name == "dihedrals")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in dihedrals table " << line << std::endl;
                    }
                }
            }
            else if(section_name == "pairs")
            {
                while(std::getline(ifs, line))
                {
                    if(line.empty() || std::regex_match(line, std::regex("^\\s*")))
                    {
                        continue; // skip empty line
                    }
                    else if(std::regex_match(line, std::regex("^;.*")))
                    {
                        continue; // skip comment line
                    }
                    else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                    {
                        next_section_name = result.str(1);
                        std::cerr << "next section name is " << next_section_name << std::endl;
                        break;
                    }
                    else
                    {
                        std::cerr << "in pairs table " << line << std::endl;
                    }
                }
            }
            else
            {
                throw std::runtime_error(
                        "[error] Unknown name table " + section_name + " detected.");
            }
        }
    }

    // read input top file
    // std::cerr << "parsing " << input_top << "..." << std::endl;

    // read input gro file
    // std::cerr << "parsing " << input_gro << "..." << std::endl;
}

Simulator read_input(int argc, char** argv)
{
    // check command line argument
    if(argc == 2)
    {
        return read_toml_input(std::string(argv[1]));
    }
    else if(argc == 4)
    {
        std::array<std::string, 3> filenames, file_suffixes;
        for(std::size_t file_idx=0; file_idx<3; ++file_idx)
        {
            filenames[file_idx]     = std::string(argv[file_idx+1]);
            file_suffixes[file_idx] = get_file_suffix(filenames[file_idx]);
        }

        std::string gro_file, top_file, itp_file;
        // get each type file name
        for(std::size_t file_idx=0; file_idx<3; ++file_idx)
        {
            if(file_suffixes[file_idx] == ".gro"){ gro_file = filenames[file_idx]; }
            if(file_suffixes[file_idx] == ".itp"){ itp_file = filenames[file_idx]; }
            if(file_suffixes[file_idx] == ".top"){ top_file = filenames[file_idx]; }
        }

        if(gro_file.empty())
        {
            throw std::runtime_error(
                    "[error] There is no `.gro` file in input files."
                    " Genesis mode needs `.gro`, `.top` and `.itp` file for input.");
        }
        if(itp_file.empty())
        {
            throw std::runtime_error(
                    "[error] There is no `.itp` file in input files."
                    " Genesis mode needs `.gro`, `.top` and `.itp` file for input.");
        }
        if(top_file.empty())
        {
            throw std::runtime_error(
                    "[error] There is no `.top` file in input files."
                    " Genesis mode needs `.gro`, `.top` and `.itp` file for input.");
        }

        read_genesis_input(gro_file, itp_file, top_file);

        throw std::runtime_error(
                "[error] Input files are " + std::string(argv[1]) + ", " +
                std::string(argv[2]) + " and " + argv[3] +
                ". Genesis mode has not yet been implemented.");
    }
    else
    {
        throw std::runtime_error(
                "Usage: " + std::string(argv[0]) + " <input.toml> or " +
                 std::string(argv[0]) + " <input.gro> <input.itp> <input.top>");
    }
}
#endif // OPEN_AICG2_PLUS_READ_INPUT_HPP
