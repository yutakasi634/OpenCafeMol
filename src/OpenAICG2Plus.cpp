#include <OpenMM.h>

#include <toml11/toml.hpp>

#include "src/util/Utility.hpp"
#include "src/util/Macro.hpp"
#include "src/input/ReadTOMLInput.hpp"
#include "src/input/ReadGenesisInput.hpp"
#include "Simulator.hpp"

#include <iostream>
#include <string>

#include <cmath>

Simulator read_input(int argc, char** argv)
{
    // check command line argument
    if(argc == 2)
    {
        const std::string file_suffix = Utility::get_file_suffix(std::string(argv[1]));
        if(file_suffix == ".toml")
        {
            return read_toml_input(std::string(argv[1]));
        }
        else if(file_suffix == ".inp")
        {
            return read_genesis_input(std::string(argv[1]));
        }
    }

    throw std::runtime_error(
            "Usage: " + std::string(argv[0]) + " <input.toml> or " +
             std::string(argv[0]) + " <input.inp>");
}

void simulate(Simulator& simulator)
{
    // excute simulation
    const auto start = std::chrono::system_clock::now();
    std::cerr << "calculation start!" << std::endl;

    simulator.run();

    const auto stop  = std::chrono::system_clock::now();
    const auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cerr << "elapsed time : ";
    std::cerr << std::fixed << std::setprecision(1);

    if(total < 1000.0)
    {
        std::cerr << total << " [msec]";
    }
    else if(total < 1000.0 * 60.0)
    {
        std::cerr << total * 0.001 << " [sec]";
    }
    else if(total < 1000.0 * 60.0 * 60.0)
    {
        std::cerr << total * 0.001 * 0.0167 << " [min]";
    }
    else if(total < 1000.0 * 60.0 * 60.0 * 24.0)
    {
        std::cerr << total * 0.001 * 0.0167 * 0.0167 << " [hr]";
    }
    else
    {
        std::cerr << total * 0.001 * 0.0167 * 0.0167 * 0.0417 << " [day]";
    }
    std::cerr << std::endl;
}


int main(int argc, char** argv)
{
    // dump library information
    std::cerr << "OpenMM Library Information" << std::endl;
    std::cerr << "    version                   : "
        + OpenMM::Platform::getOpenMMVersion() << std::endl;
    std::cerr << "    CUDA platform plugin path : "
        << OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR) << std::endl;

    // Load any shared libraries containing GPU implementations
    OpenMM::Platform::loadPluginsFromDirectory(
            OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    try {
        Simulator simulator(read_input(argc, argv));
        simulate(simulator);
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("%s\n", e.what());
        return 1; // failure!
    }
}
