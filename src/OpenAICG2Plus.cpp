#include <OpenMM.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <optional>
#include <utility>
#include "toml11/toml.hpp"
#include "Constants.hpp"
#include "Utility.hpp"

void simulateSH3(const std::string& input_file_name)
{
    // Load any shared libraries containing GPU implementations.
    OpenMM::Platform::loadPluginsFromDirectory(
        OpenMM::Platform::getDefaultPluginsDirectory());

    const std::size_t file_path_len   = input_file_name.rfind("/")+1;
    const std::size_t file_prefix_len = input_file_name.rfind(".") - file_path_len;
    const std::string file_path       = input_file_name.substr(0, file_path_len);
    const std::string file_prefix     = input_file_name.substr(file_path_len, file_prefix_len);
    // read input toml file
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

    std::cerr << "initializing forcefields..." << std::endl;
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
                std::cerr << "    BondLength    : Harmonic" << std::endl;
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
                std::cerr << "    BondLength    : Gaussian" << std::endl;
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
                std::cerr << "    BondLength    : GoContact" << std::endl;

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
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                std::cerr << "    BondAngle     : FlexibleLocalAngle" << std::endl;

                for(const AAType aa_type : AAType())
                {
                    const std::array<double, 10>& spline_table = fla_spline_table.at(aa_type);
                    OpenMM::Continuous1DFunction spline_func =
                        OpenMM::Continuous1DFunction(
                            std::vector<double>(spline_table.begin(), spline_table.end()),
                            fla_spline_min_theta, fla_spline_max_theta);
                    const double min_theta_y = spline_table[0];
                    const double max_theta_y = spline_table[9];

                    OpenMM::CustomCompoundBondForce angle_ff =
                        OpenMM::CustomCompoundBondForce(3,
                            "k * ("
                               "spline(theta)"
                               "- (step(min_theta-theta) * 30 * (theta-min_theta) + min_theta_y)"
                               "+ (step(theta-max_theta) * 30 * (theta-max_theta) + max_theta_y)"
                            ")");
                    angle_ff.addTabulatedFunction("spline", &spline_func);
                    angle_ff.addPerBondParameter("k");
                    angle_ff.addPerBondParameter("min_theta");
                    angle_ff.addPerBondParameter("min_theta_y");
                    angle_ff.addPerBondParameter("max_theta");
                    angle_ff.addPerBondParameter("max_theta_y");

                    const auto& params = toml::find<toml::array>(local_ff, "parameteres");
                    for(const auto& param : params)
                    {
                        const auto& indices =
                            toml::find<std::vector<int>>(param, "indices");
                        const double k =
                            toml::find<double>(param, "k") * OpenMM::KJPerKcal; // KJ/mol
                        angle_ff.addBond(indices,
                                {k, fla_spline_min_theta, min_theta_y,
                                 fla_spline_max_theta, max_theta_y});
                    }
                    system.addForce(&angle_ff);
                }
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                std::cerr << "    DihedralAngle : Gaussian" << std::endl;

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
                std::cerr << "    Global        : ExcludedVolume" << std::endl;
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
        simulateSH3(std::string(argv[1]));
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}
