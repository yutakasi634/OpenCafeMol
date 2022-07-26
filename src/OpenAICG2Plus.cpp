#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <optional>
#include <utility>
#include <memory>
#include "toml11/toml.hpp"
#include "src/util/Utility.hpp"
#include "src/util/Macro.hpp"
#include "input/ReadInput.hpp"
#include "Simulator.hpp"
#include "Observer.hpp"

void simulate(const std::string& input_file_name)
{
    Simulator simulator(read_input(input_file_name));

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
    // check command line argument
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.toml>" << std::endl;
        return 1;
    }

    // dump library information
    std::cerr << "OpenMM Library Information" << std::endl;
    std::cerr << "    version                   : "
        + OpenMM::Platform::getOpenMMVersion() << std::endl;
    std::cerr << "    CUDA platform plugin path : "
        << OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR) << std::endl;

    // Load any shared libraries containing GPU implementations
    OpenMM::Platform::loadPluginsFromDirectory(
            OPENAICG2PLUS_EXPAND_OPTION_STR(OPENMM_PLUGIN_DIR));

    // check CUDA platform existance
    bool cuda_platform_found = false;
    for(std::size_t idx=0; idx < OpenMM::Platform::getNumPlatforms(); ++idx)
    {
        if(OpenMM::Platform::getPlatform(idx).getName() == "CUDA")
        {
            cuda_platform_found = true;
        }
    }
    if(!cuda_platform_found)
    {
        throw std::runtime_error(
            "[error] There is no CUDA platform loaded. "
            "You need to set the correct OpenMM plugins directory path "
            "to the CMake option -DOPENMM_PLUGIN_DIR.");
    }


    try {
        simulate(std::string(argv[1]));
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}
