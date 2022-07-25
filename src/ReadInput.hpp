#ifndef OPEN_AICG2_PLUS_READ_INPUT_HPP
#define OPEN_AICG2_PLUS_READ_INPUT_HPP

#include <memory>
#include <OpenMM.h>
#include "Simulator.hpp"
#include "forcefield/HarmonicBondForceFieldGenerator.hpp"
#include "forcefield/GaussianBondForceFieldGenerator.hpp"
#include "forcefield/GoContactForceFieldGenerator.hpp"
#include "forcefield/FlexibleLocalAngleForceFieldGenerator.hpp"
#include "forcefield/GaussianDihedralForceFieldGenerator.hpp"
#include "forcefield/FlexibleLocalDihedralForceFieldGenerator.hpp"
#include "forcefield/ExcludedVolumeForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::System> read_system(const toml::value& data)
{
    auto system_ptr = std::make_unique<OpenMM::System>();

    // read systems tables
    const auto& systems   = toml::find(data, "systems");
    const auto& particles = toml::find<toml::array>(systems[0], "particles");

    std::size_t system_size = particles.size();
    std::cerr << "generating system with size " << system_size << "..." << std::endl;
    for(std::size_t i=0; i<system_size; ++i)
    {
        // set mass
        const auto& p = particles.at(i);
        system_ptr->addParticle(toml::get<double>(find_either(p, "m", "mass"))); // amu
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
                system_ptr->addForce(ff_gen.generate().release());
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
                system_ptr->addForce(ff_gen.generate().release());
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
                system_ptr->addForce(ff_gen.generate().release());
                contacted_pairs = indices_vec;
            }
            else if(interaction == "BondAngle" && potential == "FlexibleLocalAngle")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");

                for(const auto& [aa_type, spline_table] : Constant::fla_spline_table)
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
                    system_ptr->addForce(ff_gen.generate().release());
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
                system_ptr->addForce(ff_gen.generate().release());
            }
            else if(interaction == "DihedralAngle" && potential == "FlexibleLocalDihedral")
            {
                const auto& params = toml::find<toml::array>(local_ff, "parameters");

                for(const auto& [aa_pair_type, fourier_table] : Constant::fld_fourier_table)
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
                    system_ptr->addForce(ff_gen.generate().release());
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
                system_ptr->addForce(ff_gen.generate().release());
            }
        }
    }

    return system_ptr;
}

std::vector<OpenMM::Vec3> read_initial_conf(const toml::value& data)
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
        vec = toml::get<std::array<double, 3>>(find_either(p, "pos", "position"));
        initPosInNm[i] = OpenMM::Vec3(vec[0]*OpenMM::NmPerAngstrom,
                                      vec[1]*OpenMM::NmPerAngstrom,
                                      vec[2]*OpenMM::NmPerAngstrom); // location, nm
    }

    return initPosInNm;
}

Simulator read_input(const std::string input_file_name)
{
    const std::size_t file_path_len   = input_file_name.rfind("/")+1;
    const std::size_t file_prefix_len = input_file_name.rfind(".") - file_path_len;
    const std::string file_path       = input_file_name.substr(0, file_path_len);
    const std::string file_prefix     = input_file_name.substr(file_path_len, file_prefix_len);

    // read input toml file
    std::cerr << "parsing " << input_file_name << "..." << std::endl;
    auto data = toml::parse(input_file_name);

    // read system table
    const auto&    systems     = toml::find(data, "systems");
    const auto&    attr        = toml::find(systems[0], "attributes");
    const auto&    temperature = toml::find<double>(attr, "temperature");
    const std::vector<OpenMM::Vec3> initial_position_in_nm(read_initial_conf(data));

    // read simulator table
    const auto&       simulator_table = toml::find(data, "simulator");
    const std::size_t total_step      = toml::find<std::size_t>(simulator_table, "total_step");
    const std::size_t save_step       = toml::find<std::size_t>(simulator_table, "save_step");
    const double      delta_t         = toml::find<double>(simulator_table, "delta_t");

    // setup OpenMM simulator
    // In OpenMM, we cannot use different friction coefficiet, gamma, for different molecules.
    // However, cafemol use different gamma depends on the mass of each particle, and the
    // product of mass and gamma is constant in there.
    // So in this implementation, we fix the gamma to 0.2 ps^-1 temporary, correspond to
    // approximatry 0.01 in cafemol friction coefficient. We need to implement new
    // LangevinIntegrator which can use different gamma for different particles.
    return Simulator(std::move(read_system(data)),
               OpenMM::LangevinIntegrator(temperature,
                                          0.3/*friction coef ps^-1*/,
                                          delta_t*Constant::cafetime),
               initial_position_in_nm, total_step, save_step, Observer(file_prefix));
}

#endif // OPEN_AICG2_PLUS_READ_INPUT_HPP
