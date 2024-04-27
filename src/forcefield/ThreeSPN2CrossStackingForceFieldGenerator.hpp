#ifndef OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include "ThreeSPN2DefaultParameters.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

class ThreeSPN2CrossStackingForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 3>;
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    ThreeSPN2CrossStackingForceFieldGenerator(
        const ThreeSPN2CrossStackingPotentialDefaultParameter& para,
        const std::vector<indices_type>& donor_indices_vec,
        const std::vector<indices_type>& acceptor_indices_vec,
        const std::vector<std::string>&  base_kind_acceptors_vec,
        const std::pair<std::string, std::string>& base_pair,
        const index_pairs_type& ignore_list,
        const bool use_periodic)
        : donor_indices_vec_(donor_indices_vec),
          acceptor_indices_vec_(acceptor_indices_vec),
          base_kind_acceptors_vec_(base_kind_acceptors_vec), base_pair_(base_pair),
          ignore_list_(ignore_list),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("TSPN2_CS{}", ffid.gen())),
          name_      (para.name()      ),
          cutoff_    (para.cutoff()    ),
          alpha_CS_  (para.alpha_CS()  ),
          K_CS_      (para.K_CS()      ),
          K_BP_      (para.K_BP()      ),
          theta3_0_  (para.theta3_0()  ),
          epsilon_CS_(para.epsilon_CS()),
          r0_CS_     (para.r0_CS()     ),
          theta_CS_0_(para.theta_CS_0())
    {
        if(!(acceptor_indices_vec.size() == base_kind_acceptors_vec.size()))
        {
            std::ostringstream oss;
            oss << "[error] ThreeSPN2CrossStackingForceFieldGenerator: "
                   "parameter number of "
                   "acceptor_indices_vec (" << acceptor_indices_vec.size() << "), "
                   "base_kind_acceptors_vec (" << base_kind_acceptors_vec.size() << ")"
                   " is not matched "
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept override
    {
        return this->name_ + "CrossStacking "
               "(" + base_pair_.first + "-" + base_pair_.second + ")";
    }

  private:
    std::vector<indices_type>           donor_indices_vec_;
    std::vector<indices_type>           acceptor_indices_vec_;
    std::vector<std::string>            base_kind_acceptors_vec_;
    std::pair<std::string, std::string> base_pair_;
    index_pairs_type                    ignore_list_;
    bool                                use_periodic_;
    std::string                         ffgen_id_;

    std::string                   name_      ;
    double                        cutoff_    ;
    double                        alpha_CS_  ;
    double                        K_CS_      ;
    double                        K_BP_      ;
    std::map<std::string, double> theta3_0_  ;
    std::map<std::string, double> epsilon_CS_;
    std::map<std::string, double> r0_CS_     ;
    std::map<std::string, double> theta_CS_0_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_CROSS_STACKING_FORCE_FIELD_GENERATOR_HPP
