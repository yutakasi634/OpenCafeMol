#ifndef OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP

#include <regex>
#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"

std::map<std::string, std::vector<std::string>> read_genesis_itp(const std::string& input_itp)
{
    std::ifstream ifs(input_itp);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + input_itp);
    }

    std::map<std::string, std::vector<std::string>> itp_data;
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
                itp_data.insert(std::make_pair("moleculetype", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("moleculetype");
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
                        table_contents.push_back(line);
                    }
                }
            }
            else if(section_name == "atoms")
            {
                itp_data.insert(std::make_pair("atoms", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("atoms");
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
                        table_contents.push_back(line);
                    }
                }
            }
            else if(section_name == "bonds")
            {
                itp_data.insert(std::make_pair("bonds", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("bonds");
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
                        table_contents.push_back(line);
                    }
                }
            }
            else if(section_name == "angles")
            {
                itp_data.insert(std::make_pair("angles", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("angles");
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
                        table_contents.push_back(line);
                    }
                }
            }
            else if(section_name == "dihedrals")
            {
                itp_data.insert(std::make_pair("dihedrals", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("dihedrals");
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
                        table_contents.push_back(line);
                    }
                }
            }
            else if(section_name == "pairs")
            {
                itp_data.insert(std::make_pair("pairs", std::vector<std::string>{}));
                std::vector<std::string>& table_contents = itp_data.at("pairs");
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
                        table_contents.push_back(line);
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

    ifs.close();

    return itp_data;
}

std::vector<OpenMM::Vec3> read_genesis_initial_conf(const std::string& input_gro)
{
    std::ifstream ifs(input_gro);
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

std::unique_ptr<OpenMM::System> read_genesis_system(
        const std::map<std::string, std::vector<std::string>>& itp_data)
{
    auto system_ptr = std::make_unique<OpenMM::System>();

    for(auto& atoms_line : itp_data.at("atoms"))
    {
        system_ptr->addParticle(std::stof(atoms_line.substr(49, 9))); // amu
    }

    std::cerr << "generating forcefields..." << std::endl;
    // read forcefields info
    // harmonic bond force field
    std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
    std::vector<double>                              v0s;
    std::vector<double>                              ks;

    if(itp_data.count("bonds") == 1)
    {
        for(auto& bonds_line : itp_data.at("bonds"))
        {
            const std::size_t idx_i = std::stoi(bonds_line.substr( 0, 10));
            const std::size_t idx_j = std::stoi(bonds_line.substr(10, 10));
            const double      v0    = std::stof(bonds_line.substr(25, 42));
            const double      k     = std::stof(bonds_line.substr(43, 60)) * 0.5; // KJ/(mol nm^2)
            indices_vec.push_back(std::make_pair(idx_i-1, idx_j-1));
            v0s        .push_back(v0);
            ks         .push_back(k);
        }
        const auto ff_gen = HarmonicBondForceFieldGenerator(indices_vec, v0s, ks);

        std::cerr << "    BondLength    : Harmonic" << std::endl;
        system_ptr->addForce(ff_gen.generate().release());
    }

    return system_ptr;
}

Simulator read_genesis_input(
        const std::string& input_gro, const std::string& input_itp, const std::string& input_top)
{
    std::cerr << "gro file is " << input_gro << std::endl;
    std::cerr << "itp file is " << input_itp << std::endl;
    std::cerr << "top file is " << input_top << std::endl;

    // read input itp file
    std::cerr << "parsing " << input_itp << "..." << std::endl;
    std::map<std::string, std::vector<std::string>> itp_data = read_genesis_itp(input_itp);

    // read input top file
    // std::cerr << "parsing " << input_top << "..." << std::endl;

    // read input gro file
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_genesis_initial_conf(input_gro));

    return Simulator(std::move(read_genesis_system(itp_data)),
                   OpenMM::LangevinIntegrator(300.0,
                                              0.3/*friction coef ps^-1*/,
                                              0.01/* delta t*/),
                   initial_position_in_nm, 1000000, 1000, Observer("test_output"));
}

#endif // OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
