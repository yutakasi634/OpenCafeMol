#include <OpenMM.h>

#include <toml.hpp>

#include "src/util/Utility.hpp"
#include "src/util/Macro.hpp"
#include "src/util/Logger.hpp"
#include "src/input/ReadTOMLInput.hpp"
#include "src/input/ReadGenesisInput.hpp"
#include "Simulator.hpp"

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
    log_fatal("Usage: {} <input.toml> or <input.inp>", std::string(argv[0]));
}

void simulate(Simulator& simulator)
{
    // excute simulation
    const auto start = std::chrono::system_clock::now();
    log_info("calculation start!");

    simulator.run();

    const auto stop  = std::chrono::system_clock::now();
    const auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    if(total < 1000.0)
    {
        log_info("elapsed time : {} [msec]", total);
    }
    else if(total < 1000.0 * 60.0)
    {
        log_info("elapsed time : {:.1f} [sec]", total * 0.001);
    }
    else if(total < 1000.0 * 60.0 * 60.0)
    {
        log_info("elapsed time : {:.1f} [min]", total * 0.001 * 0.0167);
    }
    else if(total < 1000.0 * 60.0 * 60.0 * 24.0)
    {
        log_info("elapsed time : {:.1f} [hr]", total * 0.001 * 0.0167 * 0.0167);
    }
    else
    {
        log_info("elapsed time : {:.1f} [day]", total * 0.001 * 0.0167 * 0.0167 * 0.0417);
    }
}


int main(int argc, char** argv)
{
    // dump library information
    log_info("OpenMM Library Information");
    log_info("    version                   : {}", OpenMM::Platform::getOpenMMVersion());
    log_info("    CUDA platform plugin path : {}", OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    // Load any shared libraries containing GPU implementations
    log_info("loading OpenMM plugins...");
    OpenMM::Platform::loadPluginsFromDirectory(
            OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    for(const auto& error : OpenMM::Platform::getPluginLoadFailures())
    {
        log_warn(error);
    }

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
