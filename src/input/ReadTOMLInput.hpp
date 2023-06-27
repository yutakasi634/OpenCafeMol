#ifndef OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP

#include <memory>
#include <OpenMM.h>
#include "src/Simulator.hpp"
#include "src/SystemGenerator.hpp"
#include "src/Topology.hpp"
#include "ReadTOMLForceFieldGenerator.hpp"
#include "src/forcefield/MonteCarloAnisotropicBarostatGenerator.hpp"
#include "src/forcefield/MonteCarloMembraneBarostatGenerator.hpp"

SystemGenerator read_toml_system(const toml::value& data)
{
    // read systems table
    const auto& systems   = toml::find(data, "systems");

    // read particles info
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    Topology        topology(system_size);
    std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    std::vector<double>                     mass_vec(system_size);
    for(std::size_t idx=0; idx<system_size; ++idx)
    {
        // set mass
        const auto& p = particles.at(idx);
        mass_vec[idx] = toml::get<double>(Utility::find_either(p, "m", "mass")); // amu

        if(p.contains("group"))
        {
            const std::string group_name = toml::find<std::string>(p, "group");
            group_vec[idx] = group_name;
        }
    }
    SystemGenerator system_gen(mass_vec);

    // read boundary condition
    const auto&       simulator_table = toml::find(data, "simulator");
    const std::string boundary_type   = toml::find<std::string>(simulator_table, "boundary_type");
    bool use_periodic = false;
    if(boundary_type == "Periodic" || boundary_type == "PeriodicCuboid")
    {
        use_periodic = true;
        const auto& boundary_shape = toml::find(systems[0], "boundary_shape");
        const auto lower_bound = toml::find<std::array<double, 3>>(boundary_shape, "lower");
        const auto upper_bound = toml::find<std::array<double, 3>>(boundary_shape, "upper");
        if(lower_bound[0] != 0.0 || lower_bound[1] != 0.0 || lower_bound[2] != 0.0)
        {
            std::cerr << "[warning] Lower bound of periodic boundary box is not (0.0, 0.0, 0.0). "
                         "OpenAICG2+ only support the case lower bound is origin, "
                         "so the simulation box will be moved to satisfy this condition. "
                         "In this case, specified lower bound is ("
                      << std::fixed << std::setprecision(2)
                      << std::setw(7) << lower_bound[0] << ", "
                      << std::setw(7) << lower_bound[1] << ", "
                      << std::setw(7) << lower_bound[2] << ")" << std::endl;
        }

        const double xlength = (upper_bound[0] - lower_bound[0]) * OpenMM::NmPerAngstrom;
        const double ylength = (upper_bound[1] - lower_bound[1]) * OpenMM::NmPerAngstrom;
        const double zlength = (upper_bound[2] - lower_bound[2]) * OpenMM::NmPerAngstrom;
        system_gen.set_pbc(xlength, ylength, zlength);
    }

    // read ensemble condition
    const auto& attr = toml::find(systems[0], "attributes");
    if(attr.contains("ensemble"))
    {
        const auto& ensemble = toml::find(attr, "ensemble");
        const auto& type     = toml::find<std::string>(ensemble, "type");
        if(type == "NPT")
        {
            const auto& default_pressure =
                toml::find<std::array<double, 3>>(ensemble, "pressure");
            const auto& scale_axis =
                toml::find<std::array<bool, 3>>(ensemble, "scale_axis");
            const auto frequency = toml::find_or<std::size_t>(ensemble, "frequency", 25);

            if(!attr.contains("temperature"))
            {
                throw std::runtime_error(
                        "[error] attributes table must contains temperature in NPT ensemble case.");
            }
            const double temperature = toml::find<double>(attr, "temperature");

            auto barostat_gen =
                MonteCarloAnisotropicBarostatGenerator(
                        scale_axis, temperature, default_pressure, frequency);
            system_gen.set_barostat(
                    std::make_unique<MonteCarloAnisotropicBarostatGenerator>(barostat_gen));
        }
        else if(type == "NPgammaT" || type == "NPγT")
        {
            const auto& default_pressure = toml::find<double>(ensemble, "pressure");
            const auto& default_surface_tension =
                toml::find<double>(ensemble, "surface_tension");
            const auto frequency = toml::find_or<std::size_t>(ensemble, "frequency", 25);
            if(!attr.contains("temperature"))
            {
                throw std::runtime_error(
                        "[error] attributes table must contains temperature in NPγT"
                        " ensemble case.");
            }
            const double temperature = toml::find<double>(attr, "temperature");

            // XYMode parameter
            MonteCarloMembraneBarostatGenerator::xymode_type xymode;
            const std::string xymode_str = toml::find<std::string>(ensemble, "xymode");
            if(xymode_str == "Isotropic")
            {
                xymode = MonteCarloMembraneBarostatGenerator::xymode_type::XYIsotropic;
            }
            else if(xymode_str == "Anisotropic")
            {
                xymode = MonteCarloMembraneBarostatGenerator::xymode_type::XYAnisotropic;
            }
            else
            {
                throw std::runtime_error(
                     "[error] unknown xymode "+xymode_str+" was specified.");
            }

            // ZMode parameter
            MonteCarloMembraneBarostatGenerator::zmode_type zmode;
            const std::string zmode_str = toml::find<std::string>(ensemble, "zmode");
            if(zmode_str == "Free")
            {
                zmode = MonteCarloMembraneBarostatGenerator::zmode_type::ZFree;
            }
            else if(zmode_str == "Fixed")
            {
                zmode = MonteCarloMembraneBarostatGenerator::zmode_type::ZFixed;
            }
            else if(zmode_str == "ConstantVolume")
            {
                zmode = MonteCarloMembraneBarostatGenerator::zmode_type::ConstantVolume;
            }
            else
            {
                throw std::runtime_error(
                     "[error] unknown zmode "+xymode_str+" was specified.");
            }

            auto barostat_gen =
                MonteCarloMembraneBarostatGenerator(
                        default_pressure, default_surface_tension,
                        temperature, xymode, zmode, frequency);
            system_gen.set_barostat(
                    std::make_unique<MonteCarloMembraneBarostatGenerator>(barostat_gen));
        }
        else
        {
            throw std::runtime_error(
                    "[error] unknown ensemble type "+type+" was specified.");
        }
    }

    std::cerr << "reading forcefield tables..." << std::endl;
    // read forcefields info
    const auto ff = toml::find(data, "forcefields").at(0);
    // to distinguish same name parameters between forcefields, we need to label forcefield
    // id to that parameter. so we need to count the number of forcefield generator.
    std::size_t ffgen_count = 0;
    if(ff.contains("local"))
    {
        const auto& locals = toml::find(ff, "local").as_array();

        for(const auto& local_ff : locals)
        {
            const std::string interaction = toml::find<std::string>(local_ff, "interaction");
            const std::string potential   = toml::find<std::string>(local_ff, "potential");
            if(interaction == "BondLength" && potential == "Harmonic")
            {
                HarmonicBondForceFieldGenerator ff_gen =
                    read_toml_harmonic_bond_ff_generator(local_ff, topology, use_periodic);
                system_gen.add_ff_generator(
                        std::make_unique<HarmonicBondForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            else if(interaction == "BondLength" && potential == "Gaussian")
            {
                GaussianBondForceFieldGenerator ff_gen =
                    read_toml_gaussian_bond_ff_generator(
                        local_ff, topology, use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<GaussianBondForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            else if(interaction == "BondLength" && potential == "GoContact")
            {
                GoContactForceFieldGenerator ff_gen =
                    read_toml_go_contact_ff_generator(
                        local_ff, topology, use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<GoContactForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            else if(interaction == "BondAngle" && potential == "Harmonic")
            {
                HarmonicAngleForceFieldGenerator ff_gen =
                    read_toml_harmonic_angle_ff_generator(local_ff, topology, use_periodic);
                system_gen.add_ff_generator(
                        std::make_unique<HarmonicAngleForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                for(const auto& [aa_type, spline_table] : Constant::fla_spline_table)
                {
                    FlexibleLocalAngleForceFieldGenerator ff_gen =
                        read_toml_flexible_local_angle_ff_generator(
                            local_ff, aa_type, spline_table, topology, use_periodic, ffgen_count);
                    system_gen.add_ff_generator(
                            std::make_unique<FlexibleLocalAngleForceFieldGenerator>(ff_gen));
                }
                ++ffgen_count;
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                GaussianDihedralForceFieldGenerator ff_gen =
                    read_toml_gaussian_dihedral_ff_generator(
                            local_ff, topology, use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<GaussianDihedralForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            else if(interaction == "DihedralAngle" && potential == "FlexibleLocalDihedral")
            {
                for(const auto& [aa_pair_type, fourier_table] : Constant::fld_fourier_table)
                {
                    FlexibleLocalDihedralForceFieldGenerator ff_gen =
                        read_toml_flexible_local_dihedral_ff_generator(
                            local_ff, aa_pair_type, fourier_table, topology,
                            use_periodic, ffgen_count);
                    system_gen.add_ff_generator(
                            std::make_unique<FlexibleLocalDihedralForceFieldGenerator>(
                                ff_gen));
                    ++ffgen_count;
                }
            }
        }
    }
    topology.make_molecule("bond");

    if(ff.contains("global"))
    {
        const auto& globals = toml::find(ff, "global").as_array();

        for(const auto& global_ff : globals)
        {
            const std::string potential = toml::find<std::string>(global_ff, "potential");
            if(potential == "ExcludedVolume")
            {
                ExcludedVolumeForceFieldGenerator ff_gen =
                    read_toml_excluded_volume_ff_generator(
                        global_ff, system_size, topology, group_vec, use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<ExcludedVolumeForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            if(potential == "WCA")
            {
                if(global_ff.contains("table"))
                {
                    std::vector<std::pair<std::string, std::string>> treated_pair;
                    const auto& table = toml::find(global_ff, "table");
                    for(const auto& [first_key, first_table] : table.as_table())
                    {
                        for(const auto& [second_key, second_table] : first_table.as_table())
                        {
                            const auto name_pair = std::make_pair(first_key, second_key);
                            if(Utility::contains(treated_pair, name_pair))
                            {
                                throw std::runtime_error(
                                    "[error] parameter table for " + name_pair.first + "-" +
                                    name_pair.second + " is difined in dupulicate.");
                            }
                            else
                            {
                                const double sigma   =
                                    toml::get<double>(
                                            Utility::find_either(second_table, "sigma", "σ")) *
                                    OpenMM::NmPerAngstrom; // nm
                                const double epsilon =
                                    toml::get<double>(
                                            Utility::find_either(second_table, "epsilon", "ε")) *
                                    OpenMM::KJPerKcal; // KJPermol
                                UniformWeeksChandlerAndersenForceFieldGenerator ff_gen =
                                    read_toml_uniform_weeks_chandler_andersen_ff_generator(
                                        global_ff, system_size, sigma, epsilon, name_pair, topology,
                                        group_vec, use_periodic, ffgen_count);
                                if(ff_gen.former_group_size() == 0 || ff_gen.latter_group_size() == 0)
                                {
                                    std::cerr << "        "
                                        << "[warning] this force field generation will be skipped"
                                        << std::endl;
                                    continue;
                                }
                                system_gen.add_ff_generator(
                                    std::make_unique<
                                        UniformWeeksChandlerAndersenForceFieldGenerator>(
                                                ff_gen));
                                ++ffgen_count;
                                treated_pair.push_back(name_pair);
                            }
                        }
                    }
                }
                else
                {
                    WeeksChandlerAndersenForceFieldGenerator ff_gen =
                        read_toml_weeks_chandler_andersen_ff_generator(
                            global_ff, system_size, topology, group_vec, use_periodic, ffgen_count);
                    system_gen.add_ff_generator(
                            std::make_unique<WeeksChandlerAndersenForceFieldGenerator>(
                                ff_gen));
                    ++ffgen_count;
                }
            }
            if(potential == "DebyeHuckel")
            {
                const auto attr = toml::find(systems[0], "attributes");

                const auto temperature = toml::expect<double>(attr, "temperature");
                if(!temperature.is_ok())
                {
                    std::cerr << temperature.unwrap_err() << std::endl;
                }

                const auto ionic_strength = toml::expect<double>(attr, "ionic_strength");
                if(!ionic_strength.is_ok())
                {
                    std::cerr << ionic_strength.unwrap_err() << std::endl;
                }

                DebyeHuckelForceFieldGenerator ff_gen =
                    read_toml_debye_huckel_ff_generator(global_ff, system_size,
                        ionic_strength.unwrap(), temperature.unwrap(), topology, group_vec,
                        use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<DebyeHuckelForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            if(potential == "iSoLFAttractive")
            {
                iSoLFAttractiveForceFieldGenerator ff_gen =
                    read_toml_isolf_attractive_ff_generator(
                        global_ff, system_size, topology, group_vec, use_periodic, ffgen_count);
                system_gen.add_ff_generator(
                        std::make_unique<iSoLFAttractiveForceFieldGenerator>(ff_gen));
                ++ffgen_count;
            }
            if(potential == "LennardJonesAttractive")
            {
                if(global_ff.contains("table"))
                {
                    std::vector<std::pair<std::string, std::string>> treated_pair;
                    const auto& table = toml::find(global_ff, "table");
                    for(const auto& [first_key, first_table] : table.as_table())
                    {
                        for(const auto& [second_key, second_table] : first_table.as_table())
                        {
                            const auto name_pair = std::make_pair(first_key, second_key);
                            if(Utility::contains(treated_pair, name_pair))
                            {
                                throw std::runtime_error(
                                    "[error] parameter table for " + name_pair.first + "-" +
                                    name_pair.second + " is difined in dupulicate.");
                            }
                            else
                            {
                                const double sigma   =
                                    toml::get<double>(
                                            Utility::find_either(second_table, "sigma", "σ")) *
                                    OpenMM::NmPerAngstrom; // nm
                                const double epsilon =
                                    toml::get<double>(
                                            Utility::find_either(second_table, "epsilon", "ε")) *
                                    OpenMM::KJPerKcal; // KJ/mol
                                UniformLennardJonesAttractiveForceFieldGenerator ff_gen =
                                    read_toml_uniform_lennard_jones_attractive_ff_generator(
                                        global_ff, system_size, sigma, epsilon, name_pair,
                                        topology, group_vec, use_periodic, ffgen_count);
                                if(ff_gen.former_group_size() == 0 ||ff_gen.latter_group_size() == 0)                                {
                                    std::cerr << "        "
                                        << "[warning] this force field generation will be skipped"
                                        << std::endl;
                                    continue;
                                }
                                system_gen.add_ff_generator(
                                    std::make_unique<
                                        UniformLennardJonesAttractiveForceFieldGenerator>(
                                                ff_gen));
                                ++ffgen_count;
                                treated_pair.push_back(name_pair);
                            }
                        }
                    }
                }
                else
                {
                    throw std::runtime_error(
                        "[error] LennardJonesAttractive potential only support "
                        "table mode for providing pair-paramters.");
                }
            }
        }
    }

    if(ff.contains("external"))
    {
        const auto& externals = toml::find(ff, "external").as_array();

        for(const auto& external_ff : externals)
        {
            const std::string potential = toml::find<std::string>(external_ff, "potential");
            if(potential == "Harmonic")
            {
                const auto& env =
                    external_ff.contains("env") ? external_ff.at("env") : toml::value{};
                const auto& params = toml::find<toml::array>(external_ff, "parameters");
                for(const auto& param : params)
                {
                    HarmonicCoMPullingForceFieldGenerator ff_gen =
                        read_toml_harmonic_com_pulling_ff_generator(
                                param, use_periodic, env, ffgen_count);
                    system_gen.add_ff_generator(
                            std::make_unique<HarmonicCoMPullingForceFieldGenerator>(ff_gen));
                    ++ffgen_count;
                }
            }
        }
    }

    return system_gen;
}

std::vector<OpenMM::Vec3> read_toml_initial_conf(const toml::value& data)
{
    const auto& systems     = toml::find(data, "systems");
    const auto& particles   = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::vector<OpenMM::Vec3> initPosInNm(system_size);
    for (std::size_t i=0; i<system_size; ++i)
    {
        // set position
        const auto& p = particles.at(i);
        std::array<double, 3> vec = {0.0, 0.0, 0.0};
        vec = toml::get<std::array<double, 3>>(Utility::find_either(p, "pos", "position"));
        initPosInNm[i] = OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                      vec[1]*OpenMM::NmPerAngstrom,
                                      vec[2]*OpenMM::NmPerAngstrom); // location, nm
    }

    return initPosInNm;
}

Simulator read_toml_input(const std::string& toml_file_name)
{
    std::size_t file_path_len  = toml_file_name.rfind("/")+1;
    if(file_path_len == std::string::npos)
    {
        file_path_len = 0;
    }

    const std::string file_suffix = Utility::get_file_suffix(toml_file_name);
    if(file_suffix != ".toml")
    {
            throw std::runtime_error(
                    "[error] File suffix is `" + file_suffix + "`."
                    " Toml mode needs `.toml` file for toml.");
    }

    // read toml toml file
    std::cerr << "parsing " << toml_file_name << "..." << std::endl;
    auto data = toml::parse(toml_file_name);

    // read files table
    const auto&        files         = toml::find(data, "files");
    const auto&        output        = toml::find(files, "output");
    const std::string& output_prefix = toml::find<std::string>(output, "prefix");
    const std::string& output_path   = toml::find<std::string>(output, "path");
    const std::string& output_format = toml::find<std::string>(output, "format");
    const bool         dump_progress_bar =
        toml::find_or<bool>(output, "progress_bar", true);

    // read simulator table
    const auto&        simulator_table = toml::find(data, "simulator");
    const std::string& boundary_type   = toml::find<std::string>(simulator_table, "boundary_type");
    const std::size_t  total_step      = toml::find<std::size_t>(simulator_table, "total_step");
    const std::size_t  save_step       = toml::find<std::size_t>(simulator_table, "save_step");
    const double       delta_t         = toml::find<double>(simulator_table, "delta_t");
    std::size_t seed = 0;
    if(simulator_table.contains("seed"))
    {
        seed = toml::find<std::size_t>(simulator_table, "seed");
    }

    // read system table
    const auto& systems     = toml::find(data, "systems");
    const auto& attr        = toml::find(systems[0], "attributes");
    const auto  temperature = toml::find<double>(attr, "temperature");
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_toml_initial_conf(data));
    SystemGenerator system_gen = read_toml_system(data);

    // setup OpenMM integrator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for
    // different molecules. However, cafemol use different gamma depends on the
    // mass of each particle, and the product of mass and gamma is constant in
    // there. So in this implementation, we fix the gamma to 0.2 ps^-1 temporary,
    // correspond to approximatry 0.01 in cafemol friction coefficient. We need
    // to implement new LangevinIntegrator which can use different gamma for
    // different particles.
    auto integrator = OpenMM::LangevinIntegrator(temperature,
                                                 0.2/*friction coef ps^-1*/,
                                                 delta_t*Constant::cafetime);
    integrator.setRandomNumberSeed(seed);

    // construct observers
    const bool use_periodic =
        (boundary_type == "Periodic" || boundary_type == "PeriodicCuboid");
    std::vector<std::unique_ptr<ObserverBase>> observers;
    if(output_format == "pdb")
    {
        // read name infomation
        const std::size_t system_size = initial_position_in_nm.size();
        const auto& particles = toml::find<toml::array>(systems[0], "particles");
        std::vector<std::optional<std::string>> name_vec(system_size, std::nullopt);
        for(std::size_t idx=0; idx<system_size; ++idx)
        {
            const auto& p = particles.at(idx);
            if(p.contains("name"))
            {
                name_vec[idx] = toml::find<std::string>(p, "name");
            }
        }

        observers.push_back(
                std::make_unique<PDBObserver>(
                    output_path+output_prefix, total_step, name_vec, use_periodic));
    }
    else if(output_format == "dcd")
    {
        observers.push_back(
                std::make_unique<DCDObserver>(
                    output_path+output_prefix, total_step, save_step, delta_t, use_periodic));
    }
    else
    {
        throw std::runtime_error(
                "[error] output file format `" + output_format + "` is not supported.");
    }
    observers.push_back(std::make_unique<EnergyObserver>(output_path+output_prefix, system_gen));

    return Simulator(system_gen, integrator,
               initial_position_in_nm, total_step, save_step,
               observers, dump_progress_bar);
}

#endif // OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
