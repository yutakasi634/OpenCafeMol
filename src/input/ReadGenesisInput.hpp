#ifndef OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP

#include "src/Simulator.hpp"
#include "src/IntegratorGenerator.hpp"

#include <OpenMM.h>

std::map<std::string, std::map<std::string, std::string>> read_inp_file(const std::string& inp_file_name);

std::vector<std::string> preprocess_top_file(
        const std::string& top_file_name, const std::string& file_path);

std::map<std::string, std::vector<std::string>>
read_top_file(const std::string& topfile_name, const std::string& file_path);

std::unique_ptr<IntegratorGeneratorBase>
read_genesis_integrator_gen(
        const std::map<std::string, std::string>& dynamics,
        const std::map<std::string, std::string>& ensemble
    );

std::vector<OpenMM::Vec3> read_genesis_initial_conf(
        const std::string& input_gro, const std::string& file_path);

Simulator make_simulator_from_genesis_inputs(
        const std::map<std::string, std::map<std::string, std::string>>& inpfile_data,
        const std::string& topfile_name, const std::string& grofile_name,
        const std::string& file_path);

Simulator read_genesis_input(const std::string& inp_file_name);

#endif // OPEN_AICG2_PLUS_READ_GENESIS_INPUT_HPP
