#include <OpenMM.h>
#include <iostream>
#include "toml11/toml.hpp"

// Forward declaration of routine for printing one frame of the
// trajectory, defined later in this source file.
void writePdbFrame(int frameNum, const OpenMM::State&);
const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2);

void simulateSH3()
{
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // read input toml file
    auto data = toml::parse("input/sh3_AICG2+.toml");

    // read simulator table
    const auto& simulator         = toml::find(data, "simulator");
    const std::size_t total_step  = toml::find<std::size_t>(simulator, "total_step");
    const std::size_t save_step   = toml::find<std::size_t>(simulator, "save_step");
    const std::size_t total_frame = std::floor(total_step/save_step);

    // read systems tables
    const auto& systems     = toml::find(data, "systems");
    const auto& particles   = toml::find<toml::array>(systems[0], "particles");
    const auto& attr        = toml::find(systems[0], "attributes");
    const auto& temperature = toml::find<float>(attr, "temperature");

    OpenMM::System system;

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initPosInA(system_size);
    for (std::size_t i=0; i<system_size; ++i){
        const auto& p = particles.at(i);

        system.addParticle(toml::get<float>(find_either(p, "m", "mass"))); // mass g/mol
    }


    // setup OpenMM simulator
    OpenMM::LangevinIntegrator integrator(temperature, 0.5/*temp*/, 0.1/*temp*/);

    OpenMM::Context context(system, integrator,
            OpenMM::Platform::getPlatformByName("CudaPlatform"));

    // Set starting positions of the atoms.
    context.setPositions(initPosInA);

    // Simulate
    for (std::size_t frame_num=0; frame_num<total_frame; ++frame_num)
    {
        OpenMM::State state = context.getState(OpenMM::State::Positions);
        writePdbFrame(frame_num, state); // output coordinates

        integrator.step(save_step);
    }
}

int main()
{
    try {
        simulateSH3();
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}

// Handy homebrew PDB writer for quick-and-dirty trajectory output.
void writePdbFrame(int frameNum, const OpenMM::State& state)
{
    // Reference atomic positions in the OpenMM State.
    const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

    // Use PDB MODEL cards to number trajectory frames
    printf("MODEL     %d\n", frameNum); // start of frame
    for (int a = 0; a < (int)posInNm.size(); ++a)
    {
        printf("ATOM  %5d  AR   AR     1    ", a+1); // atom number
        printf("%8.3f%8.3f%8.3f  1.00  0.00\n",      // coordinates
            // "*10" converts nanometers to Angstroms
            posInNm[a][0]*10, posInNm[a][1]*10, posInNm[a][2]*10);
    }
    printf("ENDMDL\n"); // end of frame
}

const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2)
{
    // A functor to find a value that corresponds to either of the key.
    // If both key exists, throw an error.
    if(v.contains(key1) && v.contains(key2) != 0)
    {
        std::cerr << toml::format_error("[error] key duplicates.", v.at(key1), "here", v.at(key2),
                                        "this conflicts with the above value definition")
                 << std::endl;
    }
    if(v.contains(key1)) {return v.at(key1);}
    if(v.contains(key2)) {return v.at(key2);}

    std::cerr << toml::format_error("both keys, \"" + key1 + "\" and \"" + key2 +
                                    "\", are not found.", v, "in this table")
              << std::endl;
}
