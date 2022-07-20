#ifndef OPEN_AICG2_PLUS_OBSERVER_HPP
#define OPEN_AICG2_PLUS_OBSERVER_HPP

#include <string>
#include <iostream>
#include <OpenMM.h>

class Observer
{
  public:
    Observer(const std::string& file_prefix)
        : file_prefix_(file_prefix),
          pos_filename_("output/"+file_prefix+".pdb"),
          ene_filename_("output/"+file_prefix+".ene")
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

    void output(const std::size_t step, const OpenMM::Context& context) const
    {
        // output position
        {
            std::ofstream ofs(pos_filename_, std::ios::app);
            OpenMM::State pos = context.getState(OpenMM::State::Positions);
            writePdbFrame(ofs, step, pos);
            ofs.close();
        }

        // output energy
        {
            std::ofstream ofs(ene_filename_, std::ios::app);
            OpenMM::State ene = context.getState(OpenMM::State::Energy);
            writeEnergy(ofs, step, ene);
            ofs.close();
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

  private:
    std::string file_prefix_;
    std::string pos_filename_;
    std::string ene_filename_;
};

#endif // OPEN_AICG2_PLUS_OBSERVER_HPP
