#ifndef OPEN_AICG2_PLUS_HARMONIC_COM_PULLING_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_HARMONIC_COM_PULLING_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <limits>
#include <set>
#include <optional>
#include <OpenMM.h>
#include <fmt/core.h>
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

class HarmonicCoMPullingForceFieldGenerator : public ForceFieldGeneratorBase
{
  public:
    HarmonicCoMPullingForceFieldGenerator(
        const double k, const double v0,
        const std::vector<int>& first_group, const std::vector<int>& second_group,
        const bool use_periodic)
        : k_(k), v0_(v0), first_group_(first_group), second_group_(second_group),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("HCP{}", ffid.gen()))
    {}

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula = fmt::format(
            "0.5 * {id}_k * (distance(g1, g2) - {id}_v0)^2",
            fmt::arg("id", ffgen_id_));

        auto com_ff = std::make_unique<OpenMM::CustomCentroidBondForce>(2, potential_formula);

        com_ff->addGlobalParameter(fmt::format("{}_k", ffgen_id_), k_);
        com_ff->addGlobalParameter(fmt::format("{}_v0", ffgen_id_), v0_);
        com_ff->addGroup(first_group_);
        com_ff->addGroup(second_group_);
        com_ff->addBond({0, 1});

        com_ff->setUsesPeriodicBoundaryConditions(use_periodic_);

        return com_ff;
    }

    std::string name() const noexcept { return "HamonicCoMPulling"; }

  private:
    double           k_;
    double           v0_;
    std::vector<int> first_group_;
    std::vector<int> second_group_;
    bool             use_periodic_;
    std::string      ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_HARMONIC_COM_PULLING_FORCE_FIELD_GENERATOR_HPP
