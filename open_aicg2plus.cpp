#include <OpenMM.h>
#include <iostream>
#include <cmath>
#include "toml11/toml.hpp"

// Forward declaration of routine for printing one frame of the
// trajectory, defined later in this source file.
void writePdbFrame(int frameNum, const OpenMM::State&);
const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2);

void simulateSH3()
{
    // define physical constnt
    const double Na = 6.02214076e23;
    const double kB = 1.380649e-23; // J/K
    const double pi = 3.1415926535897932385;

    // prepare fundamental constant for calculation
    const double cafetime = std::sqrt(1.0/OpenMM::KJPerKcal)*0.1; // Ps

    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    // read input toml file
    auto data = toml::parse("input/sh3_AICG2+.toml");

    // read simulator table
    const auto& simulator              = toml::find(data, "simulator");
    const std::size_t total_step       = toml::find<std::size_t>(simulator, "total_step");
    const std::size_t save_step        = toml::find<std::size_t>(simulator, "save_step");
    const double      delta_t          = toml::find<double>(simulator, "delta_t");
    const auto&       integrator_info  = toml::find(simulator, "integrator");
    const auto&       gammas           = toml::find<toml::array>(integrator_info, "gammas");

    const std::size_t total_frame = std::floor(total_step/save_step);

    // read systems tables
    const auto& systems     = toml::find(data, "systems");
    const auto& particles   = toml::find<toml::array>(systems[0], "particles");
    const auto& attr        = toml::find(systems[0], "attributes");
    const auto& temperature = toml::find<double>(attr, "temperature");

    OpenMM::System system;

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initPosInNm(system_size);
    for (std::size_t i=0; i<system_size; ++i)
    {
        // set mass
        const auto& p = particles.at(i);
        system.addParticle(toml::get<double>(find_either(p, "m", "mass"))); // amu

        // set position
        std::array<double, 3> vec = {0.0, 0.0, 0.0};
        vec = toml::get<std::array<double, 3>>(find_either(p, "pos", "position"));
        initPosInNm[i] = OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                      vec[1]*OpenMM::NmPerAngstrom,
                                      vec[2]*OpenMM::NmPerAngstrom); // location, nm
    }

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.01 temporary. We need to implement
    // new LangevinIntegrator which can use different gamma for different particles.
    OpenMM::LangevinIntegrator integrator(temperature, 0.01/*friction coef*/, delta_t*cafetime);

    OpenMM::Context context(system, integrator,
            OpenMM::Platform::getPlatformByName("CUDA"));

    // Set starting positions of the atoms.
    context.setPositions(initPosInNm);

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
            posInNm[a][0]*OpenMM::AngstromsPerNm,
            posInNm[a][1]*OpenMM::AngstromsPerNm,
            posInNm[a][2]*OpenMM::AngstromsPerNm);
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
    if (v.contains(key1)) { return v.at(key1); }
    if(v.contains(key2))  { return v.at(key2); }

    throw std::runtime_error(toml::format_error("both keys, \"" + key1 + "\" and \"" + key2 +
                             "\", are not found.", v, "in this table"));
}

