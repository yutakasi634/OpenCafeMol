#ifndef OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_LOCAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_LOCAL_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include "ThreeSPN2DefaultParameters.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

class ThreeSPN2BasePairLocalForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 4>;

  public:
    ThreeSPN2BasePairLocalForceFieldGenerator(
        const ThreeSPN2BasePairPotentialParameterBase& para,
        const std::vector<indices_type>& indices_vec,
        const std::pair<std::string, std::string>& base_pair,
        const bool use_periodic)
        : indices_vec_(indices_vec),
          base_pair_(base_pair),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("TSPN2BPL{}", ffid.gen())),
          alpha_BP_  (para.alpha_BP  ()),
          K_BP_      (para.K_BP      ()),
          epsilon_BP_(para.epsilon_BP()),
          r0_        (para.r0        ()),
          theta0_1_  (para.theta0_1  ()),
          theta0_2_  (para.theta0_2  ()),
          phi0_      (para.phi0      ()),
          name_      (para.name      ())
    {}

    std::unique_ptr<OpenMM::Force> generate() const override;

    std::string name() const override
    {
        return name_+"BasePairLocal "
               "(" + base_pair_.first + "-" + base_pair_.second + ")";
    }

  private:
    std::vector<indices_type> indices_vec_;
    std::pair<std::string, std::string> base_pair_;
    bool                      use_periodic_;
    std::string               ffgen_id_;

    double                        alpha_BP_;
    double                        K_BP_;
    std::map<std::string, double> epsilon_BP_;
    std::map<std::string, double> r0_;
    std::map<std::string, double> theta0_1_;
    std::map<std::string, double> theta0_2_;
    std::map<std::string, double> phi0_;
    std::string                   name_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_LOCAL_FORCE_FIELD_GENERATOR_HPP
