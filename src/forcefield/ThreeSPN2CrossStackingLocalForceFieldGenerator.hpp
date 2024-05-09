#ifndef OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_LOCAL_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_LOCAL_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include "ThreeSPN2DefaultParameters.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

class ThreeSPN2CrossStackingLocalForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 6>;

  public:
    ThreeSPN2CrossStackingLocalForceFieldGenerator(
        const ThreeSPN2CrossStackingPotentialParameterBase& para,
        const std::vector<indices_type>& indices_vec,
        const std::vector<std::string>&  base_kind_vec,
        const std::pair<std::string, std::string>& base_pair,
        const bool use_periodic)
        : indices_vec_(indices_vec),
          base_kind_vec_(base_kind_vec), base_pair_(base_pair),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("TSPN2_CSL{}", ffid.gen())),
          name_      (para.name()      ),
          alpha_CS_  (para.alpha_CS()  ),
          K_CS_      (para.K_CS()      ),
          K_BP_      (para.K_BP()      ),
          theta3_0_  (para.theta3_0()  ),
          epsilon_CS_(para.epsilon_CS()),
          r0_CS_     (para.r0_CS()     ),
          theta_CS_0_(para.theta_CS_0())

    {
        if(!(indices_vec.size() == base_kind_vec.size()))
        {
            std::ostringstream oss;
            oss << "[error] ThreeSPN2CrossStackingLocalForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec ("   << indices_vec.size() << "), "
                   "base_kind_vec (" << base_kind_vec.size() << ")"
                   " is not matched "
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept override
    {
        return this->name_ + "CrossStackingLocal "
               "(" + base_pair_.first + "-" + base_pair_.second + ")";
    }

  private:
    std::vector<indices_type>           indices_vec_;
    std::vector<std::string>            base_kind_vec_;
    std::pair<std::string, std::string> base_pair_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;

    std::string                   name_      ;
    double                        alpha_CS_  ;
    double                        K_CS_      ;
    double                        K_BP_      ;
    std::map<std::string, double> theta3_0_  ;
    std::map<std::string, double> epsilon_CS_;
    std::map<std::string, double> r0_CS_     ;
    std::map<std::string, double> theta_CS_0_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_LOCAL_FORCE_FIELD_GENERATOR_HPP
