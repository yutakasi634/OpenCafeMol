#ifndef OPEN_AICG2_PLUS_OBSERVER_HPP
#define OPEN_AICG2_PLUS_OBSERVER_HPP

#include <string>
#include <iostream>
#include <OpenMM.h>
#include "util/ProgressBar.hpp"

class Observer
{
  public:
    Observer(const std::string& file_prefix, bool output_progress = true)
        : pos_filename_("output/"+file_prefix+".pdb"),
          ene_filename_("output/"+file_prefix+".ene"),
          output_progress_(output_progress),
          progress_bar_(/* width of bar = */ 50)
    {
        std::cerr << "initializing observer..." << std::endl;

        this->clear_file(pos_filename_);
        std::cerr << "    output trajectory file : " << pos_filename_ << std::endl;


        this->clear_file(ene_filename_);
        std::cerr << "    output energy file     : " << ene_filename_ << std::endl;
        {
            std::ofstream ofs(ene_filename_, std::ios::out);
            ofs << "# unit of energy : kcal/mol" << std::endl;
            ofs << "# timestep potential_energy kinetic_energy" << std::endl;
            ofs.close();
        }
    }

    void initialize(const std::size_t total_step)
    {
        total_step_ = total_step;
    }

    void output(const std::size_t step, const OpenMM::Context& context) const
    {
        // output position
        {
            std::ofstream ofs(pos_filename_, std::ios::app);
            OpenMM::State pos = context.getState(OpenMM::State::Positions);
            write_pdb_frame(ofs, step, pos);
            ofs.close();
        }

        // output energy
        {
            std::ofstream ofs(ene_filename_, std::ios::app);
            OpenMM::State ene = context.getState(OpenMM::State::Energy);
            write_energy(ofs, step, ene);
            ofs.close();
        }

        if(output_progress_)
        {
            progress_bar_.format(step, total_step_, std::cerr);
        }
    }

    void finalize() const
    {
        if(output_progress_)
        {
            progress_bar_.finalize(std::cerr);
        }
    }

  private:
    void clear_file(const std::string& filename) const
    {
        std::ofstream ofs(filename);
        if(not ofs.good())
        {
            throw std::runtime_error("file open error : " + filename);
        }
        ofs.close();
        return;
    }

    // Handy homebrew PDB writer for quick-and-dirty trajectory output.
    void write_pdb_frame(std::ofstream& fp, int frame_num, const OpenMM::State& state) const
    {
        // Reference atomic positions in the OpenMM State.
        const std::vector<OpenMM::Vec3>& pos_in_nm = state.getPositions();

        // Use PDB MODEL cards to number trajectory frames
        fp << "MODEL     " << frame_num << std::endl; // start of frame
        for (int a = 0; a < (int)pos_in_nm.size(); ++a)
        {
            fp << std::setprecision(3);
            fp << "ATOM  " << std::setw(5) << a+1 << "  AR   AR     1    "; // atom number
            fp << std::setw(8) << std::fixed
               << std::setw(8) << std::fixed << pos_in_nm[a][0]*OpenMM::AngstromsPerNm
               << std::setw(8) << std::fixed << pos_in_nm[a][1]*OpenMM::AngstromsPerNm
               << std::setw(8) << std::fixed << pos_in_nm[a][2]*OpenMM::AngstromsPerNm << "  1.00  0.00" << std::endl;
        }
        fp << "ENDMDL" << std::endl; // end of frame
    }

    void write_energy(std::ofstream& fp, int frame_num, const OpenMM::State& state) const
    {
        const double pot_ene = state.getPotentialEnergy() * OpenMM::KcalPerKJ; // kcal/mol
        const double kin_ene = state.getKineticEnergy() * OpenMM::KcalPerKJ; // kcal/mol
        fp << std::setw(11) << std::left << frame_num << ' ';
        fp << std::setw(16) << std::right << std::fixed << pot_ene;
        fp << "  " << std::setw(14) << std::right << std::fixed << kin_ene << std::endl;
    }

  private:
    std::string pos_filename_;
    std::string ene_filename_;
    bool        output_progress_;
    ProgressBar progress_bar_;
    std::size_t total_step_;
};

#endif // OPEN_AICG2_PLUS_OBSERVER_HPP
