#include "SystemGenerator.hpp"
#include "util/Logger.hpp"

void SystemGenerator::add_ff_generator(std::unique_ptr<ForceFieldGeneratorBase>&& ff_gen_ptr)
{
    const std::size_t group_id = ffname_groupid_map_.size();
    const std::string ff_name  = ff_gen_ptr->name();
    if(ff_name.find("FlexibleLocalAngle") != std::string::npos)
    {
        if(ffname_groupid_map_.find("FlexibleLocalAngle") ==
                ffname_groupid_map_.end())
        {
            ffname_groupid_map_.insert(
                    std::make_pair("FlexibleLocalAngle", group_id));
        }
    }
    else if(ff_name.find("FlexibleLocalDihedral") != std::string::npos)
    {
        if(ffname_groupid_map_.find("FlexibleLocalDihedral") ==
                ffname_groupid_map_.end())
        {
            ffname_groupid_map_.insert(
                    std::make_pair("FlexibleLocalDihedral", group_id));
        }
    }
    else
    {
        ffname_groupid_map_.insert(std::make_pair(ff_name, group_id));
    }

    ff_gen_ptrs_.push_back(std::move(ff_gen_ptr));
}

std::unique_ptr<OpenMM::System> SystemGenerator::generate() const
{
    std::unique_ptr system_ptr = std::make_unique<OpenMM::System>();

    log_info("generating system with size {}...", mass_vec_.size());
    for(auto& mass : mass_vec_)
    {
        system_ptr->addParticle(mass);
    }

    // set boundary condition
    if(edge_lengthes_opt_)
    {
        const std::array<double, 3>& edge_length = edge_lengthes_opt_.value();
        log_info("    boundary type is periodic cuboid");
        log_info("        the box size is {:7.2f} x {:7.2f} x {:7.2f} nm^3",
                edge_length.at(0), edge_length.at(1), edge_length.at(2));

        system_ptr->setDefaultPeriodicBoxVectors(
                OpenMM::Vec3(edge_length[0],            0.0,            0.0),
                OpenMM::Vec3(           0.0, edge_length[1],            0.0),
                OpenMM::Vec3(           0.0,            0.0, edge_length[2]));
    }
    else
    {
        log_info("    boundary type is unlimited");
    }

    // set forcefields
    log_info("generating forcefields...");
    for(auto& ff_gen_ptr : ff_gen_ptrs_)
    {
        const std::string ff_name = ff_gen_ptr->name();
        std::unique_ptr<OpenMM::Force> ff_ptr = ff_gen_ptr->generate();

        // set forcegroup id
        if(ff_name.find("FlexibleLocalAngle") != std::string::npos)
        {
            ff_ptr->setForceGroup(ffname_groupid_map_.at("FlexibleLocalAngle"));
        }
        else if(ff_name.find("FlexibleLocalDihedral") != std::string::npos)
        {
            ff_ptr->setForceGroup(ffname_groupid_map_.at("FlexibleLocalDihedral"));
        }
        else if(ff_name.find("CombinatorialGoContact") != std::string::npos)
        {
            ff_ptr->setForceGroup(ffname_groupid_map_.at("CombinatorialGoContact"));
        }
        else
        {
            ff_ptr->setForceGroup(ffname_groupid_map_.at(ff_name));
        }
        system_ptr->addForce(ff_ptr.release());
    }

    // set barostat
    if(barostat_gen_opt_)
    {
        log_info("generating barostat...");
        const auto& barostat_gen_ptr = barostat_gen_opt_.value();

        log_assert(edge_lengthes_opt_.has_value(),
            "barostat should be used with periodic boundary condition.");

        log_info("    barostat is {}", barostat_gen_ptr->name());
        const std::size_t frequency = barostat_gen_ptr->frequency();
        log_info("        pressure change frequency : {:14}", frequency);

        const double temperature = barostat_gen_ptr->temperature();
        log_info("        barostat temperature      : {:14f} K", temperature);

        barostat_gen_ptr->dump_info();

        system_ptr->addForce(barostat_gen_ptr->generate().release());
    }

    return system_ptr;
}

