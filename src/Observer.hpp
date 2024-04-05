#ifndef OPEN_AICG2_PLUS_OBSERVER_HPP
#define OPEN_AICG2_PLUS_OBSERVER_HPP

#include <string>
#include <iostream>
#include <OpenMM.h>
#include "src/SystemGenerator.hpp"
#include "util/Utility.hpp"
#include "util/ProgressBar.hpp"

class ObserverBase
{
  public:
    virtual void initialize(const std::unique_ptr<OpenMM::System>&) = 0;
    virtual void output(const std::size_t /*step*/, const OpenMM::Context& /*context*/) = 0;
    virtual void finalize() const = 0;
    virtual std::string name() const = 0;
};

class PDBObserver final : public ObserverBase
{
  public:
    PDBObserver(const std::string& file_prefix, const std::size_t total_step,
                const std::vector<std::optional<std::string>> name_vec, const bool use_periodic)
        : pos_filename_(file_prefix+".pdb"), total_step_(total_step), use_periodic_(use_periodic)
    {
        Utility::clear_file(pos_filename_);
        std::cerr << "    output trajectory file    : " << pos_filename_ << std::endl;

        const std::size_t system_size = name_vec.size();
        name_vec_.resize(system_size, "UNK");
        for(std::size_t idx=0; idx<name_vec.size(); ++idx)
        {
            const auto& name = name_vec[idx];
            if(name)
            {
                name_vec_[idx] = name.value();
            }
        }
    }

    void initialize(const std::unique_ptr<OpenMM::System>&) override { return; }

    void output(const std::size_t step, const OpenMM::Context& context) override
    {
        // output position
        std::ofstream ofs(pos_filename_, std::ios::app);
        OpenMM::State pos = context.getState(OpenMM::State::Positions,
                                             use_periodic_);
        write_pdb_frame(ofs, step, pos);
        ofs.close();
    }

    void finalize() const override { return; }
    std::string name() const { return "PDBObserver"; }

  private:
    std::string              pos_filename_;
    std::size_t              total_step_;
    bool                     use_periodic_;
    std::vector<std::string> name_vec_;

  private:
    // Handy homebrew PDB writer for quick-and-dirty trajectory output.
    void write_pdb_frame(std::ofstream& fp, int frame_num, const OpenMM::State& state) const
    {
        // Reference atomic positions in the OpenMM State.
        const std::vector<OpenMM::Vec3>& pos_in_nm = state.getPositions();

        // Use PDB MODEL cards to number trajectory frames
        fp << "MODEL     " << frame_num << std::endl; // start of frame
        for (std::size_t idx=0; idx<pos_in_nm.size(); ++idx)
        {
            const std::string& name = name_vec_[idx];
            fp << std::setprecision(3);
            fp << "ATOM  "                         //                          1-6
               << std::setw(5) << idx+1 << " "     // atom serial number       7-12
               << "  CA "                          // atom name               13-17
               << std::setw(3) << name << "  "     // residue name            18-22
               << "   1    " ;                     // residue sequence number 23-30
            fp << std::setw(8) << std::fixed
               << std::setw(8) << std::fixed << pos_in_nm[idx][0]*OpenMM::AngstromsPerNm
               << std::setw(8) << std::fixed << pos_in_nm[idx][1]*OpenMM::AngstromsPerNm
               << std::setw(8) << std::fixed << pos_in_nm[idx][2]*OpenMM::AngstromsPerNm << "  1.00  0.00" << std::endl;
        }
        fp << "ENDMDL" << std::endl; // end of frame
    }
};



class DCDObserver final : public ObserverBase
{
  public:
    DCDObserver(const std::string& file_prefix, const std::size_t total_step,
        const std::size_t save_interval, const float delta_t, const bool use_periodic)
        : pos_filename_(file_prefix+"_position.dcd"),
          vel_filename_(file_prefix+"_velocity.dcd"),
          total_step_(total_step), delta_t_(delta_t),
          save_interval_(save_interval), use_periodic_(use_periodic)
    {
        Utility::clear_file(pos_filename_);
        std::cerr << "    output trajectory file    : " << pos_filename_ << std::endl;
        Utility::clear_file(vel_filename_);
        std::cerr << "    output velocity file      : " << vel_filename_ << std::endl;
    }

    void initialize(const std::unique_ptr<OpenMM::System>& system_ptr) override
    {
        // position file
        std::ofstream pos_ofs(pos_filename_, std::ios::binary | std::ios::app);
        write_dcd_header(pos_ofs, system_ptr);
        pos_ofs.close();

        // velocity file
        std::ofstream vel_ofs(vel_filename_, std::ios::binary | std::ios::app);
        write_dcd_header(vel_ofs, system_ptr);
        vel_ofs.close();

        // buffer to convert dcd coordinate format
        const std::size_t num_of_particles(system_ptr->getNumParticles());
        buffer_pos_x_.resize(num_of_particles);
        buffer_pos_y_.resize(num_of_particles);
        buffer_pos_z_.resize(num_of_particles);
        buffer_vel_x_.resize(num_of_particles);
        buffer_vel_y_.resize(num_of_particles);
        buffer_vel_z_.resize(num_of_particles);

        return;
    }

    void output(const std::size_t step, const OpenMM::Context& context) override
    {
        // output position
        std::ofstream pos_ofs(pos_filename_, std::ios::binary | std::ios::app);
        OpenMM::State pos = context.getState(OpenMM::State::Positions, use_periodic_);
        write_dcd_frame(pos_ofs, pos);
        pos_ofs.close();

        // output velocity
        std::ofstream vel_ofs(vel_filename_, std::ios::binary | std::ios::app);
        OpenMM::State vel = context.getState(OpenMM::State::Velocities, use_periodic_);
        write_dcd_velocity(vel_ofs, vel);
        vel_ofs.close();

        return;
    }

    void finalize() const override { return; }
    std::string name() const { return "DCDObserver"; }

  private:
    std::string        pos_filename_;
    std::string        vel_filename_;
    std::size_t        total_step_;
    float              delta_t_;
    std::size_t        save_interval_;
    bool               use_periodic_;
    std::vector<float> buffer_pos_x_;
    std::vector<float> buffer_pos_y_;
    std::vector<float> buffer_pos_z_;
    std::vector<float> buffer_vel_x_;
    std::vector<float> buffer_vel_y_;
    std::vector<float> buffer_vel_z_;

  private:
    void write_dcd_header(
            std::ofstream& ofs, const std::unique_ptr<OpenMM::System>& system_ptr) const
    {
        // write first block
        {
            const std::int32_t block_size(84);
            Utility::write_as_bytes(ofs, block_size);

            ofs.write("CORD", 4);

            const std::int32_t total_frame(std::floor(total_step_ / save_interval_));
            Utility::write_as_bytes(ofs, total_frame);

            const std::int32_t index_of_first(0);
            Utility::write_as_bytes(ofs, index_of_first);

            const std::int32_t save_interval_i32(save_interval_);
            Utility::write_as_bytes(ofs, save_interval_i32);

            const std::int32_t total_step_i32(total_step_);
            Utility::write_as_bytes(ofs, total_step_i32);

            const std::int32_t total_chains(0); // TODO: put correct molecule num
            Utility::write_as_bytes(ofs, total_chains);

            const std::int32_t zero(0);
            // 4 * integers with null flag
            for(std::size_t i=0; i<4; ++i) {Utility::write_as_bytes(ofs, zero);}

            Utility::write_as_bytes(ofs, delta_t_);

            const std::int32_t has_unitcell = use_periodic_;
            Utility::write_as_bytes(ofs, has_unitcell);

            // 8 * integers with null flag
            for(std::size_t i=0; i<8; ++i) {Utility::write_as_bytes(ofs, zero);}

            const std::int32_t version(24);
            Utility::write_as_bytes(ofs, version);

            Utility::write_as_bytes(ofs, block_size);
        }

        // write second block
        {
            const std::int32_t block_size(84);
            Utility::write_as_bytes(ofs, block_size);

            const std::int32_t number_of_lines(1);
            Utility::write_as_bytes(ofs, number_of_lines);

            const char comment[80] = "OpenMoCa -- written by Yutaka Murata 2022";

            ofs.write(comment, 80);

            Utility::write_as_bytes(ofs, block_size);
        }

        // write thrid block
        {
            const std::int32_t block_size(4);
            Utility::write_as_bytes(ofs, block_size);

            const std::int32_t number_of_particles(system_ptr->getNumParticles());
            Utility::write_as_bytes(ofs, number_of_particles);

            Utility::write_as_bytes(ofs, block_size);
        }
    }

    void write_unitcell(std::ofstream& ofs, const OpenMM::State& state)
    {
        OpenMM::Vec3 A, B, C;
        state.getPeriodicBoxVectors(A, B, C);

        // This software only use cuboid box for periodic boundary condition.
        // So, each A, B, C vector corresponds to X, Y, Z axis.
        const double x_A   = A[0] * OpenMM::AngstromsPerNm;
        const double y_B   = B[1] * OpenMM::AngstromsPerNm;
        const double z_C   = C[2] * OpenMM::AngstromsPerNm;
        const double alpha = 90.0;
        const double beta  = 90.0;
        const double gamma = 90.0;

        const std::int32_t block_size = sizeof(double) * 6;
        Utility::write_as_bytes(ofs, block_size);

        Utility::write_as_bytes(ofs, x_A);
        Utility::write_as_bytes(ofs, gamma);
        Utility::write_as_bytes(ofs, y_B);
        Utility::write_as_bytes(ofs, beta);
        Utility::write_as_bytes(ofs, alpha);
        Utility::write_as_bytes(ofs, z_C);

        Utility::write_as_bytes(ofs, block_size);
    }

    void write_dcd_frame(std::ofstream& ofs, const OpenMM::State& state)
    {
        // Reference atomic position in the OpenMM State.
        const std::vector<OpenMM::Vec3>& pos_in_nm = state.getPositions();

        assert(this->buffer_pos_x_.size() == pos_in_nm.size());
        assert(this->buffer_pos_y_.size() == pos_in_nm.size());
        assert(this->buffer_pos_z_.size() == pos_in_nm.size());

        // write periodic box information
        if(use_periodic_) { write_unitcell(ofs, state); }

        // write position
        {
            for(std::size_t idx=0; idx<pos_in_nm.size(); ++idx)
            {
                this->buffer_pos_x_[idx] = static_cast<float>(pos_in_nm[idx][0]*OpenMM::AngstromsPerNm);
                this->buffer_pos_y_[idx] = static_cast<float>(pos_in_nm[idx][1]*OpenMM::AngstromsPerNm);
                this->buffer_pos_z_[idx] = static_cast<float>(pos_in_nm[idx][2]*OpenMM::AngstromsPerNm);
            }
            const std::int32_t block_size(sizeof(float) * pos_in_nm.size());
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_pos_x_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_pos_y_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_pos_z_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
        }
        return;
    }

    void write_dcd_velocity(std::ofstream& ofs, const OpenMM::State& state)
    {
        // Reference atomic velocity in the OpenMM State.
        const std::vector<OpenMM::Vec3>& vel_in_nmps = state.getVelocities();

        assert(this->buffer_vel_x_.size() == vel_in_nmps.size());
        assert(this->buffer_vel_y_.size() == vel_in_nmps.size());
        assert(this->buffer_vel_z_.size() == vel_in_nmps.size());

        // write velocity
        {
            for(std::size_t idx=0; idx<vel_in_nmps.size(); ++idx)
            {
                this->buffer_vel_x_[idx] = static_cast<float>(vel_in_nmps[idx][0]*OpenMM::AngstromsPerNm*Constant::cafetime);
                this->buffer_vel_y_[idx] = static_cast<float>(vel_in_nmps[idx][1]*OpenMM::AngstromsPerNm*Constant::cafetime);
                this->buffer_vel_z_[idx] = static_cast<float>(vel_in_nmps[idx][2]*OpenMM::AngstromsPerNm*Constant::cafetime);
            }
            const std::int32_t block_size(sizeof(float) * vel_in_nmps.size());
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_vel_x_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_vel_y_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_vel_z_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
        }
        return;
    }
};

class EnergyObserver final : public ObserverBase
{
  public:
    EnergyObserver(const std::string& file_prefix, const SystemGenerator& system_gen)
    : ene_filename_(file_prefix+".ene"),
      ffname_groupid_map_(system_gen.ffname_groupid_map())
    {
        Utility::clear_file(ene_filename_);
        std::cerr << "    output energy file        : " << ene_filename_ << std::endl;
    }

    void initialize(const std::unique_ptr<OpenMM::System>&) override
    {
        std::ofstream ofs(ene_filename_, std::ios::out);
        ofs << "# unit of energy : kcal/mol" << std::endl;
        ofs << "# timestep";
        for(auto ffname_groupid : ffname_groupid_map_)
        {
            ofs << " " << std::setw(24) << ffname_groupid.first;
        }
        ofs << " kinetic_energy" << std::endl;
        ofs.close();
        return;
    }

    void output(const std::size_t step, const OpenMM::Context& context) override
    {
        // output energy
        std::ofstream ofs(ene_filename_, std::ios::app);
        write_energy(ofs, step, context);
        ofs.close();
    }

    void finalize() const override { return; }
    std::string name() const { return "EnergyObserver"; }

  private:
    void write_energy(std::ofstream& fp, int frame_num,
                      const OpenMM::Context& context) const
    {
        fp << std::setw(10) << std::left  << frame_num;
        for(auto ffname_groupid : ffname_groupid_map_)
        {
            const std::size_t groupid = ffname_groupid.second;
            const OpenMM::State state =
                context.getState(OpenMM::State::Energy, false, int(1)<<groupid);
            const double pot_ene = state.getPotentialEnergy() * OpenMM::KcalPerKJ; // kcal/mol
            fp << ' ' << std::setw(24) << std::right << std::fixed << pot_ene;
        }
        const OpenMM::State state = context.getState(OpenMM::State::Energy);
        const double kin_ene = state.getKineticEnergy() * OpenMM::KcalPerKJ; // kcal/mol
        fp << ' ' << std::setw(14) << std::right << std::fixed << kin_ene << std::endl;
    }

  private:
    std::string                        ene_filename_;
    std::map<std::string, std::size_t> ffname_groupid_map_;
};

#endif // OPEN_AICG2_PLUS_OBSERVER_HPP
