#ifndef OPEN_AICG2_PLUS_OBSERVER_HPP
#define OPEN_AICG2_PLUS_OBSERVER_HPP

#include <string>
#include <iostream>
#include <OpenMM.h>
#include "util/Utility.hpp"
#include "util/ProgressBar.hpp"

class ObserverBase
{
  public:
    virtual void initialize(const std::unique_ptr<OpenMM::System>&) = 0;
    virtual void output(const std::size_t /*step*/, const OpenMM::Context& /*context*/) = 0;
    virtual void finalize() const = 0;
    virtual const std::string name() const = 0;
};

class PDBObserver final : public ObserverBase
{
  public:
    PDBObserver(const std::string& file_prefix, const std::size_t total_step,
                const std::vector<std::optional<std::string>> name_vec)
        : pos_filename_(file_prefix+".pdb"), total_step_(total_step)
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
                                             /*enforcePeriodicBox*/ true);
        write_pdb_frame(ofs, step, pos);
        ofs.close();
    }

    void finalize() const override { return; }
    const std::string name() const { return "PDBObserver"; }

  private:
    std::string              pos_filename_;
    std::size_t              total_step_;
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
        const std::size_t save_interval, const float delta_t)
        : dcd_filename_(file_prefix+".dcd"), total_step_(total_step), delta_t_(delta_t),
          save_interval_(save_interval)
    {
        Utility::clear_file(dcd_filename_);
        std::cerr << "    output trajectory file    : " << dcd_filename_ << std::endl;
    }

    void initialize(const std::unique_ptr<OpenMM::System>& system_ptr) override
    {
        std::ofstream ofs(dcd_filename_, std::ios::binary | std::ios::app);
        write_dcd_header(ofs, system_ptr);
        ofs.close();

        // buffer to convert dcd coordinate format
        const std::size_t num_of_particles(system_ptr->getNumParticles());
        buffer_x_.resize(num_of_particles);
        buffer_y_.resize(num_of_particles);
        buffer_z_.resize(num_of_particles);
        return;
    }

    void output(const std::size_t step, const OpenMM::Context& context) override
    {
        // output position
        std::ofstream ofs(dcd_filename_, std::ios::binary | std::ios::app);
        OpenMM::State pos = context.getState(OpenMM::State::Positions,
                                             /*enforcePeriodicBox*/ true);
        write_dcd_frame(ofs, pos);
        ofs.close();
        return;
    }

    void finalize() const override { return; }
    const std::string name() const { return "DCDObserver"; }

  private:
    std::string        dcd_filename_;
    std::size_t        total_step_;
    float              delta_t_;
    std::size_t        save_interval_;
    std::vector<float> buffer_x_;
    std::vector<float> buffer_y_;
    std::vector<float> buffer_z_;

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

            // 9 * integers with null flag
            for(std::size_t i=0; i<9; ++i) {Utility::write_as_bytes(ofs, zero);}

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

            const char comment[80] = "OpenAICG2+ -- written by Yutaka Murata 2022";

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

        // TODO: implement 4th block existing case
        return;
    }

    void write_dcd_frame(std::ofstream& ofs, const OpenMM::State& state)
    {
        // Reference atomic position in the OpenMM State.
        const std::vector<OpenMM::Vec3>& pos_in_nm = state.getPositions();

        assert(this->buffer_x_.size() == pos_in_nm.size());
        assert(this->buffer_y_.size() == pos_in_nm.size());
        assert(this->buffer_z_.size() == pos_in_nm.size());

        // write position
        {
            for(std::size_t idx=0; idx<pos_in_nm.size(); ++idx)
            {
                this->buffer_x_[idx] = static_cast<float>(pos_in_nm[idx][0]*OpenMM::AngstromsPerNm);
                this->buffer_y_[idx] = static_cast<float>(pos_in_nm[idx][1]*OpenMM::AngstromsPerNm);
                this->buffer_z_[idx] = static_cast<float>(pos_in_nm[idx][2]*OpenMM::AngstromsPerNm);
            }
            const std::int32_t block_size(sizeof(float) * pos_in_nm.size());
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_x_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_y_.data()),
                          block_size);
                Utility::write_as_bytes(ofs, block_size);
            }
            {
                Utility::write_as_bytes(ofs, block_size);
                ofs.write(reinterpret_cast<const char*>(this->buffer_z_.data()),
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
    EnergyObserver(const std::string& file_prefix) : ene_filename_(file_prefix+".ene")
    {
        Utility::clear_file(ene_filename_);
        std::cerr << "    output energy file        : " << ene_filename_ << std::endl;
    }

    void initialize(const std::unique_ptr<OpenMM::System>&) override
    {
        std::ofstream ofs(ene_filename_, std::ios::out);
        ofs << "# unit of energy : kcal/mol" << std::endl;
        ofs << "# timestep potential_energy kinetic_energy" << std::endl;
        ofs.close();
        return;
    }

    void output(const std::size_t step, const OpenMM::Context& context) override
    {
        // output energy
        std::ofstream ofs(ene_filename_, std::ios::app);
        OpenMM::State ene = context.getState(OpenMM::State::Energy);
        write_energy(ofs, step, ene);
        ofs.close();
    }

    void finalize() const override { return; }
    const std::string name() const { return "EnergyObserver"; }

  private:
    void write_energy(std::ofstream& fp, int frame_num, const OpenMM::State& state) const
    {
        const double pot_ene = state.getPotentialEnergy() * OpenMM::KcalPerKJ; // kcal/mol
        const double kin_ene = state.getKineticEnergy() * OpenMM::KcalPerKJ; // kcal/mol
        fp << std::setw(11) << std::left << frame_num << ' ';
        fp << std::setw(16) << std::right << std::fixed << pot_ene;
        fp << "  " << std::setw(14) << std::right << std::fixed << kin_ene << std::endl;
    }

  private:
    std::string ene_filename_;
};

#endif // OPEN_AICG2_PLUS_OBSERVER_HPP
