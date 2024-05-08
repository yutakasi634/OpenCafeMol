#ifndef OPEN_AICG2_PLUS_SYSTEM_HPP
#define OPEN_AICG2_PLUS_SYSTEM_HPP

#include "forcefield/ForceFieldGeneratorBase.hpp"
#include "forcefield/BarostatGeneratorBase.hpp"

#include <OpenMM.h>

#include <array>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

class SystemGenerator
{
  public:
    SystemGenerator(const std::vector<double> mass_vec)
    : mass_vec_(mass_vec), edge_lengthes_opt_(std::nullopt), barostat_gen_opt_(std::nullopt)
    {}

    void add_ff_generator(std::unique_ptr<ForceFieldGeneratorBase>&& ff_gen_ptr)
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

    void set_barostat(std::unique_ptr<BarostatGeneratorBase>&& baro_gen)
    {
        barostat_gen_opt_ = std::move(baro_gen);
    }

    void set_pbc(const double xlength, const double ylength, const double zlength)
    {
        edge_lengthes_opt_ = {xlength, ylength, zlength}; // unit of edge length is Nm
    }

    std::unique_ptr<OpenMM::System> generate() const
    {
        std::unique_ptr system_ptr = std::make_unique<OpenMM::System>();

        std::cerr << "generating system with size " << mass_vec_.size()
                  << "..." << std::endl;
        for(auto& mass : mass_vec_)
        {
            system_ptr->addParticle(mass);
        }

        // set boundary condition
        if(edge_lengthes_opt_)
        {
            const std::array<double, 3>& edge_length = edge_lengthes_opt_.value();
            std::cerr << "    boundary type is periodic cuboid" << std::endl;
            std::cerr << "        the box size is " << std::fixed << std::setprecision(2)
                      << std::setw(7) << edge_length[0] << " Nm x "
                      << std::setw(7) << edge_length[1] << " Nm x "
                      << std::setw(7) << edge_length[2] << " Nm" << std::endl;
            system_ptr->setDefaultPeriodicBoxVectors(
                    OpenMM::Vec3(edge_length[0],            0.0,            0.0),
                    OpenMM::Vec3(           0.0, edge_length[1],            0.0),
                    OpenMM::Vec3(           0.0,            0.0, edge_length[2]));
        }
        else
        {
            std::cerr << "    boundary type is unlimited" << std::endl;
        }

        // set forcefields
        std::cerr << "generating forcefields..." << std::endl;
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
            std::cerr << "generating barostat..." << std::endl;
            const auto& barostat_gen_ptr = barostat_gen_opt_.value();
            if(!edge_lengthes_opt_)
            {
                throw std::runtime_error(
                        "[error] barostat should be used with periodic boundary condition.");
            }

            std::cerr << "    barostat is " << barostat_gen_ptr->name() << std::endl;
            const std::size_t frequency = barostat_gen_ptr->frequency();
            std::cerr << "        pressure change frequency : "
                      << std::setw(14) << frequency << std::endl;

            const double temperature = barostat_gen_ptr->temperature();
            std::cerr << "        barostat temperature      : "
                      << std::setw(14) << std::fixed << temperature  << " K" << std::endl;

            barostat_gen_ptr->dump_info();

            system_ptr->addForce(barostat_gen_ptr->generate().release());
        }

        return system_ptr;
    }

    std::map<std::string, std::size_t> ffname_groupid_map() const
    {
        return ffname_groupid_map_;
    }

  private:
    std::vector<double>                                   mass_vec_;
    std::vector<std::unique_ptr<ForceFieldGeneratorBase>> ff_gen_ptrs_;

    // for periodic boundary condition
    std::optional<std::array<double, 3>>                  edge_lengthes_opt_;

    // for barostat
    std::optional<std::unique_ptr<BarostatGeneratorBase>> barostat_gen_opt_;

    // for EnergyObserver
    std::map<std::string, std::size_t>                    ffname_groupid_map_;
};

#endif // OPEN_AICG2_PLUS_SYSTEM_HPP
