#ifndef OPEN_AICG2_PLUS_SYSTEM_HPP
#define OPEN_AICG2_PLUS_SYSTEM_HPP

#include "forcefield/ForceFieldGeneratorBase.hpp"
#include "forcefield/BarostatGeneratorBase.hpp"

#include <OpenMM.h>

#include <array>
#include <optional>
#include <string>
#include <vector>

class SystemGenerator
{
  public:

    SystemGenerator(std::vector<double> mass_vec)
        : mass_vec_(std::move(mass_vec)),
          edge_lengthes_opt_(std::nullopt),
          barostat_gen_opt_(std::nullopt)
    {}

    void add_ff_generator(std::unique_ptr<ForceFieldGeneratorBase>&& ff_gen_ptr);

    void set_barostat(std::unique_ptr<BarostatGeneratorBase>&& baro_gen)
    {
        barostat_gen_opt_ = std::move(baro_gen);
    }

    void set_pbc(const double xlength, const double ylength, const double zlength)
    {
        edge_lengthes_opt_ = {xlength, ylength, zlength}; // unit of edge length is Nm
    }

    std::unique_ptr<OpenMM::System> generate() const;

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
