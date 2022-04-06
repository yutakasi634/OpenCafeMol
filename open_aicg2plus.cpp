#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <optional>
#include <utility>
#include "toml11/toml.hpp"
#include "math/constants.hpp"

// Forward declaration of routine for printing one frame of the
// trajectory, defined later in this source file.
void writePdbFrame(std::ofstream& fp, int frameNum, const OpenMM::State&);
void writeEnergy  (std::ofstream& fp, int frameNum, const OpenMM::State&);
const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2);

void simulateSH3()
{
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    std::string file_prefix = "sh3_AICG2+";
    // read input toml file
    auto data = toml::parse("input/" + file_prefix + ".toml");

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

    // for exclusion list of Excluded Volume
    std::vector<std::pair<std::size_t, std::size_t>> exclusion_pairs;

    std::cerr << "initializing forcefields" << std::endl;
    // read forcefields info
    const auto  ff = toml::find(data, "forcefields").at(0);
    if(ff.contains("local"))
    {
        const auto& locals = toml::find(ff, "local").as_array();

        for(const auto& local_ff : locals)
        {
            const std::string interaction = toml::find<std::string>(local_ff, "interaction");
            const std::string potential   = toml::find<std::string>(local_ff, "potential");
            if(interaction == "BondLength" && potential == "Harmonic")
            {
                std::cerr << "    BondLegth interaction : Harmonic potential" << std::endl;
                OpenMM::HarmonicBondForce* bond_ff = new OpenMM::HarmonicBondForce();

                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                for(const auto& param : params)
                {
                    const auto&  indices =
                        toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
                    const double v0 =
                        toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
                    const double k =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal *
                        OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm; // KJ/(mol nm^2)
                    bond_ff->addBond(indices.first, indices.second, v0, k);

                    exclusion_pairs.push_back(std::make_pair(indices.first, indices.second));
                }
                system.addForce(bond_ff);
            }
            else if(interaction == "BondLength" && potential == "Gaussian")
            {
                std::cerr << "    BondLegth interaction : Gaussian potential" << std::endl;
                OpenMM::CustomBondForce* bond_ff =
                    new OpenMM::CustomBondForce("k*exp(-(r-v0)^2/(2*sigma^2))");
                bond_ff->addPerBondParameter("k");
                bond_ff->addPerBondParameter("v0");
                bond_ff->addPerBondParameter("sigma");

                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                for(const auto& param : params)
                {
                    const auto& indices =
                        toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
                    const double k  =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                    const double v0 =
                        toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
                    const double sigma =
                        toml::get<double>(
                                find_either(param, "sigma", "σ")) * OpenMM::NmPerAngstrom; // nm
                    bond_ff->addBond(indices.first, indices.second, {k, v0, sigma});

                    exclusion_pairs.push_back(std::make_pair(indices.first, indices.second));
                }
                system.addForce(bond_ff);
            }
            else if(interaction == "BondLength" && potential == "GoContact")
            {
                std::cerr << "    BondLegth interaction : GoContact potential" << std::endl;

                // TODO: enable to optimization based on cutoff
                OpenMM::CustomBondForce* contact_ff =
                    new OpenMM::CustomBondForce("k*(5*(r0/r)^12-6*(r0/r)^10)");
                contact_ff->addPerBondParameter("k");
                contact_ff->addPerBondParameter("r0");

                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                for(const auto& param : params)
                {
                    const auto& indices =
                        toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
                    const double k =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                    const double r0 =
                        toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
                    contact_ff->addBond(indices.first, indices.second, {k, r0});

                    exclusion_pairs.push_back(std::make_pair(indices.first, indices.second));
                }
                system.addForce(contact_ff);
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                std::cerr << "    DihedralAngle interaction : Gaussian potential" << std::endl;
                OpenMM::CustomTorsionForce* torsion_ff =
                    new OpenMM::CustomTorsionForce("k*exp(-(theta-theta0)^2/(2*sigma^2))");
                torsion_ff->addPerTorsionParameter("k");
                torsion_ff->addPerTorsionParameter("theta0");
                torsion_ff->addPerTorsionParameter("sigma");

                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                for(const auto& param : params)
                {
                    const auto& indices =
                        toml::find<std::array<std::size_t, 4>>(param, "indices");
                    const double k  =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                    const double theta0 = toml::find<double>(param, "v0"); // radiuns
                    const double sigma =
                        toml::get<double>(find_either(param, "sigma", "σ")); // radiuns
                    torsion_ff->addTorsion(
                            indices[0], indices[1], indices[2], indices[3], {k, theta0, sigma});

                    exclusion_pairs.push_back(std::make_pair(indices[0], indices[3]));
                }
                system.addForce(torsion_ff);
            }
        }
    }

    if(ff.contains("global"))
    {
        const auto& globals = toml::find(ff, "global").as_array();

        for(const auto& global_ff : globals)
        {
            const std::string potential = toml::find<std::string>(global_ff, "potential");
            if(potential == "ExcludedVolume")
            {
                std::cerr << "    ExcludedVolume potential" << std::endl;
                // TODO: add cutoff
                const std::string exv_expression = "epsilon*((sigma1+sigma2)/r)^12";
                OpenMM::CustomNonbondedForce* exv_ff =
                    new OpenMM::CustomNonbondedForce(exv_expression);
                exv_ff->addPerParticleParameter("sigma");

                const double eps =
                    toml::find<double>(global_ff, "epsilon") * OpenMM::KJPerKcal; // KJPermol
                exv_ff->addGlobalParameter("epsilon", eps);

                const auto&  params = toml::find<toml::array>(global_ff, "parameters");
                std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
                for(const auto& param : params)
                {
                    const std::size_t index  = toml::find<std::size_t>(param, "index");
                    const double      radius =
                        toml::find<double>(param, "radius") * OpenMM::NmPerAngstrom; // nm
                    radius_vec.at(index) = radius;
                }

                const auto& itr = std::find(radius_vec.begin(), radius_vec.end(), std::nullopt);
                if(itr != radius_vec.end())
                {
                    std::size_t index = std::distance(radius_vec.begin(), itr);
                    throw std::runtime_error("[error] particle index " + std::to_string(index) +
                            " does not have radius parameter.");
                }

                for(const auto& radius : radius_vec)
                {
                    exv_ff->addParticle({ radius.value() });
                }

                for(const auto& exclusion_pair : exclusion_pairs)
                {
                    exv_ff->addExclusion(exclusion_pair.first, exclusion_pair.second);
                }

                system.addForce(exv_ff);
            }
        }
    }

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.2 ps^-1 temporary, correspond to
    // approximatry 0.01 in cafemol friction coefficient. We need to implement new
    // LangevinIntegrator which can use different gamma for different particles.
    std::cerr << "initializing integrator with " << std::endl;
    std::cerr << "    temperature : "
        << std::setw(7) << std::fixed << std::setprecision(2) << temperature << " K" << std::endl;
    std::cerr << "    delta t     : "
        << std::setw(7) << std::fixed << std::setprecision(2) << delta_t << " cafetime" << std::endl;
    OpenMM::LangevinIntegrator integrator(
            temperature, 0.3/*friction coef ps^-1*/, delta_t*cafetime);

    OpenMM::Context context(system, integrator,
            OpenMM::Platform::getPlatformByName("CUDA"));

    // Set starting positions of the atoms.
    context.setPositions(initPosInNm);

    std::ofstream pdb_fp("output/" + file_prefix + ".pdb", std::ios::out);
    std::ofstream ene_fp("output/" + file_prefix + ".ene" , std::ios::out);
    ene_fp << "# unit of energy : kcal/mol" << std::endl;
    ene_fp << "# timestep  potential_energy  kinetic_energy" << std::endl;

    // Simulate
    for (std::size_t frame_num=0; frame_num<total_frame; ++frame_num)
    {
        OpenMM::State pos    = context.getState(OpenMM::State::Positions);
        OpenMM::State energy = context.getState(OpenMM::State::Energy);
        writePdbFrame(pdb_fp, frame_num*save_step, pos); // output coordinates
        writeEnergy  (ene_fp, frame_num*save_step, energy); // output energy

        integrator.step(save_step);
    }
    pdb_fp.close();
    ene_fp.close();
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
void writePdbFrame(std::ofstream& fp, int frameNum, const OpenMM::State& state)
{
    // Reference atomic positions in the OpenMM State.
    const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

    // Use PDB MODEL cards to number trajectory frames
    fp << "MODEL     " << frameNum << std::endl; // start of frame
    for (int a = 0; a < (int)posInNm.size(); ++a)
    {
        fp << std::setprecision(3);
        fp << "ATOM  " << std::setw(5) << a+1 << "  AR   AR     1    "; // atom number
        fp << std::setw(8) << std::fixed
           << std::setw(8) << std::fixed << posInNm[a][0]*OpenMM::AngstromsPerNm
           << std::setw(8) << std::fixed << posInNm[a][1]*OpenMM::AngstromsPerNm
           << std::setw(8) << std::fixed << posInNm[a][2]*OpenMM::AngstromsPerNm << "  1.00  0.00" << std::endl;
    }
    fp << "ENDMDL" << std::endl; // end of frame
}

void writeEnergy(std::ofstream& fp, int frameNum, const OpenMM::State& state)
{
    const double pot_ene = state.getPotentialEnergy() * OpenMM::KcalPerKJ; // kcal/mol
    const double kin_ene = state.getKineticEnergy() * OpenMM::KcalPerKJ; // kcal/mol
    fp << std::setw(11) << std::left << frameNum << ' ';
    fp << std::setw(16) << std::right << std::fixed << pot_ene;
    fp << "  " << std::setw(14) << std::right << std::fixed << kin_ene << std::endl;
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

