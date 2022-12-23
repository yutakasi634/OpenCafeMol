#ifndef OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP

#include <memory>
#include <OpenMM.h>
#include "src/Simulator.hpp"
#include "src/Topology.hpp"
#include "ReadTOMLForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::System> read_toml_system(const toml::value& data)
{
    auto system_ptr = std::make_unique<OpenMM::System>();

    // read systems tables
    const auto& systems   = toml::find(data, "systems");
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    Topology topology(system_size);
    std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    std::cerr << "generating system with size " << system_size << "..." << std::endl;
    for(std::size_t idx=0; idx<system_size; ++idx)
    {
        // set mass
        const auto& p = particles.at(idx);
        system_ptr->addParticle(toml::get<double>(Utility::find_either(p, "m", "mass"))); // amu
        if(p.contains("group"))
        {
            const std::string group_name = toml::find<std::string>(p, "group");
            group_vec[idx] = group_name;
        }
    }

    std::cerr << "generating forcefields..." << std::endl;
    // read forcefields info
    const auto ff = toml::find(data, "forcefields").at(0);
    if(ff.contains("local"))
    {
        const auto& locals = toml::find(ff, "local").as_array();

        for(const auto& local_ff : locals)
        {
            const std::string interaction = toml::find<std::string>(local_ff, "interaction");
            const std::string potential   = toml::find<std::string>(local_ff, "potential");
            if(interaction == "BondLength" && potential == "Harmonic")
            {
                const auto ff_gen =
                    read_toml_harmonic_bond_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondLength" && potential == "Gaussian")
            {
                const auto ff_gen =
                    read_toml_gaussian_bond_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondLength" && potential == "GoContact")
            {
                const auto ff_gen =
                    read_toml_go_contact_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondAngle" && potential == "Harmonic")
            {
                const auto ff_gen =
                    read_toml_harmonic_angle_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                for(const auto& [aa_type, spline_table] : Constant::fla_spline_table)
                {
                    const auto ff_gen =
                        read_toml_flexible_local_angle_ff_generator(
                                local_ff, aa_type, spline_table, topology);
                    system_ptr->addForce(ff_gen.generate().release());
                }
            }
            else if(interaction == "DihedralAngle" && potential == "Gaussian")
            {
                const auto ff_gen =
                    read_toml_gaussian_dihedral_ff_generator(local_ff, topology);
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "DihedralAngle" && potential == "FlexibleLocalDihedral")
            {
                for(const auto& [aa_pair_type, fourier_table] : Constant::fld_fourier_table)
                {
                    const auto ff_gen =
                        read_toml_flexible_local_dihedral_ff_generator(
                            local_ff, aa_pair_type, fourier_table, topology);
                    system_ptr->addForce(ff_gen.generate().release());
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
                const auto ff_gen =
                    read_toml_excluded_volume_ff_generator(
                            global_ff, system_size, topology, group_vec);
                system_ptr->addForce(ff_gen.generate().release());
            }
            if(potential == "WCA")
            {
                const auto ff_gen =
                    read_toml_weeks_chandler_andersen_ff_generator(
                            global_ff, system_size, topology, group_vec);
                system_ptr->addForce(ff_gen.generate().release());
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

                const auto ff_gen =
                    read_toml_debye_huckel_ff_generator(global_ff, system_size,
                        ionic_strength.unwrap(), temperature.unwrap(), topology, group_vec);
                system_ptr->addForce(ff_gen.generate().release());
            }
            if(potential == "iSoLFAttractive")
            {
                const auto ff_gen =
                    read_toml_isolf_attractive_ff_generator(
                            global_ff, system_size, topology, group_vec);
                system_ptr->addForce(ff_gen.generate().release());
            }
        }
    }

    return system_ptr;
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

    // read system table
    const auto& systems     = toml::find(data, "systems");
    const auto& attr        = toml::find(systems[0], "attributes");
    const auto& temperature = toml::find<double>(attr, "temperature");
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_toml_initial_conf(data));

    // read simulator table
    const auto&       simulator_table = toml::find(data, "simulator");
    const std::size_t total_step      = toml::find<std::size_t>(simulator_table, "total_step");
    const std::size_t save_step       = toml::find<std::size_t>(simulator_table, "save_step");
    const double      delta_t         = toml::find<double>(simulator_table, "delta_t");

    // construct observers
    std::vector<std::unique_ptr<ObserverBase>> observers;
    if(output_format == "pdb")
    {
        observers.push_back(
                std::make_unique<PDBObserver>(
                    output_path+output_prefix, total_step));
    }
    else if(output_format == "dcd")
    {
        observers.push_back(
                std::make_unique<DCDObserver>(
                    output_path+output_prefix, total_step, save_step, delta_t));
    }
    else
    {
        throw std::runtime_error(
                "[error] output file format `" + output_format + "` is not supported.");
    }
    observers.push_back(std::make_unique<EnergyObserver>(output_path+output_prefix));

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.2 ps^-1 temporary, correspond to
    // approximatry 0.01 in cafemol friction coefficient. We need to implement new
    // LangevinIntegrator which can use different gamma for different particles.
    return Simulator(std::move(read_toml_system(data)),
               OpenMM::LangevinIntegrator(temperature,
                                          0.2/*friction coef ps^-1*/,
                                          delta_t*Constant::cafetime),
               initial_position_in_nm, total_step, save_step, observers);
}

#endif // OPEN_AICG2_PLUS_READ_TOML_INPUT_HPP
