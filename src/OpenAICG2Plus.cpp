#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <optional>
#include <utility>
#include <memory>
#include "toml11/toml.hpp"
#include "util/Utility.hpp"
#include "util/Macro.hpp"
#include "forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "forcefield/GaussianBondForceFieldGenerator.hpp"
#include "forcefield/GoContactForceFieldGenerator.hpp"
#include "forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "forcefield/ExcludedVolumeForceFieldGenerator.hpp"

void simulate(const std::string& input_file_name)
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

    const std::size_t file_path_len   = input_file_name.rfind("/")+1;
    const std::size_t file_prefix_len = input_file_name.rfind(".") - file_path_len;
    const std::string file_path       = input_file_name.substr(0, file_path_len);
    const std::string file_prefix     = input_file_name.substr(file_path_len, file_prefix_len);
    // read input toml file
    std::cerr << "parsing " << input_file_name << "..." << std::endl;
    auto data = toml::parse(input_file_name);

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
    std::cerr << "generating system with size " << system_size << "..." << std::endl;
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
    std::vector<std::pair<std::size_t, std::size_t>> bonded_pairs;
    std::vector<std::pair<std::size_t, std::size_t>> contacted_pairs;

    std::cerr << "generating forcefields..." << std::endl;
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
                auto bond_ff = std::make_unique<OpenMM::HarmonicBondForce>();

                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
                std::vector<double>                              v0s;
                std::vector<double>                              ks;

                for(const auto& param : params)
                {
                    const auto&  indices =
                        toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
                    const double v0 =
                        toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm
                    const double k =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal *
                        OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm; // KJ/(mol nm^2)

                    indices_vec.push_back(indices);
                    v0s        .push_back(v0);
                    ks         .push_back(k);
                }
                const auto ff_gen = HarmonicBondForceFieldGenerator(indices_vec, v0s, ks);
                system.addForce(ff_gen.generate().release());
                bonded_pairs = ff_gen.indices();
            }
            else if(interaction == "BondLength" && potential == "Gaussian")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
                std::vector<double>                              v0s;
                std::vector<double>                              ks;
                std::vector<double>                              sigmas;

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

                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                    v0s        .push_back(v0);
                    sigmas     .push_back(sigma);
                }
                const auto ff_gen = GaussianBondForceFieldGenerator(indices_vec, ks, v0s, sigmas);
                system.addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondLength" && potential == "GoContact")
            {
                // TODO: enable to optimization based on cutoff
                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                std::vector<std::pair<std::size_t, std::size_t>> indices_vec;
                std::vector<double>                              ks;
                std::vector<double>                              r0s;

                for(const auto& param : params)
                {
                    const auto& indices =
                        toml::find<std::pair<std::size_t, std::size_t>>(param, "indices");
                    const double k =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                    const double r0 =
                        toml::find<double>(param, "v0") * OpenMM::NmPerAngstrom; // nm

                    ks         .push_back(k);
                    indices_vec.push_back(indices);
                    r0s        .push_back(r0);
                }
                const auto ff_gen = GoContactForceFieldGenerator(indices_vec, ks, r0s);
                system.addForce(ff_gen.generate().release());
                contacted_pairs = indices_vec;
            }
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");

                for(const auto& [aa_type, spline_table] : fla_spline_table)
                {
                    std::vector<std::vector<std::size_t>> indices_vec;
                    std::vector<double>                   ks;

                    for(const auto& param : params)
                    {
                        const std::string y = toml::find<std::string>(param, "y");
                        if(y.substr(3, 3) == aa_type) // y is like "y1_PHE"
                        {
                            const auto& indices =
                                toml::find<std::vector<std::size_t>>(param, "indices");
                            const double k =
                                toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol

                            indices_vec.push_back(indices);
                            ks         .push_back(k);
                        }
                    }
                    const auto ff_gen = FlexibleLocalAngleForceFieldGenerator(
                                            indices_vec, ks, spline_table, aa_type);
                    system.addForce(ff_gen.generate().release());
                }
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");
                std::vector<std::array<std::size_t, 4>> indices_vec;
                std::vector<double>                     ks;
                std::vector<double>                     theta0s;
                std::vector<double>                     sigmas;

                for(const auto& param : params)
                {
                    const auto& indices =
                        toml::find<std::array<std::size_t, 4>>(param, "indices");
                    const double k  =
                        toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                    const double theta0 = toml::find<double>(param, "v0"); // radiuns
                    const double sigma =
                        toml::get<double>(find_either(param, "sigma", "σ")); // radiuns

                    indices_vec.push_back(indices);
                    ks         .push_back(k);
                    theta0s    .push_back(theta0);
                    sigmas     .push_back(sigma);
                }
                const auto ff_gen = GaussianDihedralForceFieldGenerator(indices_vec, ks, theta0s, sigmas);
                system.addForce(ff_gen.generate().release());
            }
            else if(interaction == "DihedralAngle" && potential == "FlexibleLocalDihedral")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");

                for(const auto& [aa_pair_type, fourier_table] : fld_fourier_table)
                {
                    std::vector<std::array<std::size_t, 4>> indices_vec;
                    std::vector<double>                     ks;
                    std::string                             aa_pair_name;

                    if(aa_pair_type.second == "GLY") // R1-GLY case
                    {
                        for(const auto& param : params)
                        {
                            const std::string coef    = toml::find<std::string>(param, "coef");
                            const auto&       indices =
                                toml::find<std::array<std::size_t, 4>>(param, "indices");
                            const double      k       = toml::find<double>(param, "k");

                            if(coef.substr(4, 3) == "GLY")
                            {
                                indices_vec.push_back(indices);
                                ks         .push_back(k);
                            }
                        }
                        aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
                    }
                    else if(aa_pair_type.second == "PRO")
                    {
                        if(aa_pair_type.first == "GLY") // GLY-PRO case
                        {
                            for(const auto& param : params)
                            {
                                const std::string coef    =
                                    toml::find<std::string>(param, "coef");
                                const auto&       indices =
                                    toml::find<std::array<std::size_t, 4>>(param, "indices");
                                const double      k       = toml::find<double>(param, "k");

                                if(coef == "GLY-PRO")
                                {
                                    indices_vec.push_back(indices);
                                    ks         .push_back(k);
                                }
                            }
                            aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
                        }
                        else
                        {
                            for(const auto& param : params)
                            {
                                const std::string coef    =
                                    toml::find<std::string>(param, "coef");
                                const auto&       indices =
                                    toml::find<std::array<std::size_t, 4>>(param, "indices");
                                const double      k       = toml::find<double>(param, "k");

                                if(coef.substr(4, 3) == "PRO") // R2-PRO case
                                {
                                    indices_vec.push_back(indices);
                                    ks         .push_back(k);
                                }
                            }
                            aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
                        }
                    }
                    else // R1-R3 case
                    {
                        for(const auto& param : params)
                        {
                            const std::string coef    =
                                toml::find<std::string>(param, "coef");
                            const auto&       indices =
                                toml::find<std::array<std::size_t, 4>>(param, "indices");
                            const double      k       = toml::find<double>(param, "k");

                            const std::string second_aa = coef.substr(4, 3);
                            if(second_aa != "GLY" && second_aa != "PRO")
                            {
                                if(coef.substr(0, 3) == aa_pair_type.first)
                                {
                                    indices_vec.push_back(indices);
                                    ks         .push_back(k);
                                }
                            }
                        }
                        aa_pair_name = aa_pair_type.first + "-" + aa_pair_type.second;
                    }
                    const auto ff_gen = FlexibleLocalDihedralForceFieldGenerator(
                            indices_vec, ks, fourier_table, aa_pair_name);
                    system.addForce(ff_gen.generate().release());
                }
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
                const double eps =
                    toml::find<double>(global_ff, "epsilon") * OpenMM::KJPerKcal; // KJPermol
                const double cutoff = toml::find_or(global_ff, "cutoff", 2.0);

                const auto&  params = toml::find<toml::array>(global_ff, "parameters");
                std::vector<std::optional<double>> radius_vec(system_size, std::nullopt);
                for(const auto& param : params)
                {
                    const std::size_t index  = toml::find<std::size_t>(param, "index");
                    const double      radius =
                        toml::find<double>(param, "radius") * OpenMM::NmPerAngstrom; // nm
                    radius_vec.at(index) = radius;
                }

                const auto ff_gen = ExcludedVolumeForceFieldGenerator(
                        eps, cutoff, radius_vec, bonded_pairs, contacted_pairs);
                system.addForce(ff_gen.generate().release());
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
    std::cerr << "initializing integrator..." << std::endl;
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

    std::cerr << "output file information" << std::endl;
    std::string output_coordinate_file = "output/" + file_prefix + ".pdb";
    std::cerr << "    output trajectory file : " << output_coordinate_file << std::endl;
    std::ofstream pdb_fp(output_coordinate_file, std::ios::out);
    if(not pdb_fp.good())
    {
        throw std::runtime_error("file open error : " + output_coordinate_file);
    }

    std::string output_energy_file     = "output/" + file_prefix + ".ene";
    std::cerr << "    output energy file     : " << output_energy_file << std::endl;
    std::ofstream ene_fp(output_energy_file, std::ios::out);
    if(not ene_fp.good())
    {
        throw std::runtime_error("file open error : " + output_coordinate_file);
    }

    ene_fp << "# unit of energy : kcal/mol" << std::endl;
    ene_fp << "# timestep  potential_energy  kinetic_energy" << std::endl;

    // Simulate
    const auto start = std::chrono::system_clock::now();
    std::cerr << "calculation start!" << std::endl;

    for (std::size_t frame_num=0; frame_num<total_frame; ++frame_num)
    {
        OpenMM::State pos    = context.getState(OpenMM::State::Positions);
        OpenMM::State energy = context.getState(OpenMM::State::Energy);
        writePdbFrame(pdb_fp, frame_num*save_step, pos); // output coordinates
        writeEnergy  (ene_fp, frame_num*save_step, energy); // output energy

        integrator.step(save_step);
    }

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

    pdb_fp.close();
    ene_fp.close();
}

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.toml>" << std::endl;
        return 1;
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
