#ifndef OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP

#include "src/Simulator.hpp"
#include "src/SystemGenerator.hpp"
#include "src/IntegratorGenerator.hpp"

#include <toml_fwd.hpp>
#include <OpenMM.h>

#include <memory>

SystemGenerator read_toml_system(const toml::value& data);

std::unique_ptr<IntegratorGeneratorBase>
read_toml_integrator_gen(const toml::value& root);

std::vector<OpenMM::Vec3> read_toml_initial_conf(const toml::value& data);
std::vector<OpenMM::Vec3> read_toml_initial_vel(const toml::value& data);

Simulator read_toml_input(const std::string& toml_file_name);

#endif // OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
