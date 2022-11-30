#ifndef OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP

#include <regex>
#include "src/forcefield/HarmonicBondForceFieldGenerator.hpp"

std::map<std::string, std::map<std::string, std::string>> read_inp_file(const std::string& inp_file_name)
{
    std::cerr << "reading " << inp_file_name << "..." << std::endl;
    std::ifstream ifs(inp_file_name);
    if(!ifs.good())
    {
        throw std::runtime_error("[error] file open error ->" + inp_file_name);
    }

    std::map<std::string, std::map<std::string, std::string>> inp_data;
    std::string line, section_name, next_section_name;
    while(std::getline(ifs, line))
    {
        std::smatch result;
        if(!next_section_name.empty())
        {
            section_name = next_section_name;
            inp_data.insert(std::make_pair(section_name, std::map<std::string, std::string>()));
            next_section_name.clear();
        }
        else if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
        {
            section_name = result.str(1);
            inp_data.insert(std::make_pair(section_name, std::map<std::string, std::string>()));
        }

        if(!section_name.empty())
        {
            std::map<std::string, std::string>& table_contents = inp_data.at(section_name);
            while(std::getline(ifs, line))
            {
                if(std::regex_match(line, result, std::regex("\\s*\\[\\s*(\\S.*?)\\s*]\\s*")))
                {
                    next_section_name = result.str(1);
                    break;
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
        const std::string& topfile_name, const std::string& grofile_name, const std::string& file_path)
{
    auto system_ptr = std::make_unique<OpenMM::System>();

    const std::map<std::string, std::vector<std::string>>& top_data =
        read_top_file(topfile_name, file_path);

    if(top_data.find("atoms") == top_data.end())
    {
        throw std::runtime_error(
                "[error] There is no atoms section. Top file must contain one atoms section.");
    }
    for(auto& atoms_line : top_data.at("atoms"))
    {
        system_ptr->addParticle(std::stof(atoms_line.substr(49, 9))); // amu
    }

    std::cerr << "generating forcefields..." << std::endl;
    if(top_data.find("bonds") != top_data.end())
    {
        std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
        std::vector<double>                              v0s;
        std::vector<double>                              ks;

        for(auto& bonds_line : top_data.at("bonds"))
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

    const std::vector<OpenMM::Vec3> initial_position_in_nm(
            read_genesis_initial_conf(grofile_name, file_path));

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.2 ps^-1 temporary, correspond to
    // approximatry 0.01 in cafemol friction coefficient. We need to implement new
    // LangevinIntegrator which can use different gamma for different particles.
    return Simulator(std::move(system_ptr),
                   OpenMM::LangevinIntegrator(300.0,
                                              0.3/*friction coef ps^-1*/,
                                              0.01/* delta t*/),
                   initial_position_in_nm, 1000000, 1000, Observer("test_output"));
}

Simulator read_inp_input(const std::string& inp_file_name)
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

    return make_simulator_from_genesis_inputs(topfile_name, grofile_name, file_path);
}

#endif // OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
