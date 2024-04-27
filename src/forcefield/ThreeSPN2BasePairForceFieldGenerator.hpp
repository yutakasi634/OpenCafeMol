#ifndef OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include "ThreeSPN2DefaultParameters.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

class ThreeSPN2BasePairForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type     = std::array<std::size_t, 2>;
    using index_pairs_type = std::vector<std::pair<std::size_t, std::size_t>>;

  public:
    ThreeSPN2BasePairForceFieldGenerator(
        const ThreeSPN2BasePairPotentialDefaultParameter& para,
        const std::vector<indices_type>& donor_indices_vec,
        const std::vector<indices_type>& acceptor_indices_vec,
        const std::pair<std::string, std::string>& base_pair,
        const index_pairs_type& ignore_list,
        const bool use_periodic)
        : donor_indices_vec_(donor_indices_vec),
          acceptor_indices_vec_(acceptor_indices_vec),
          base_pair_(base_pair), ignore_list_(ignore_list),
          use_periodic_(use_periodic),
          ffgen_id_(fmt::format("TSPN2BP{}", ffid.gen())),
          cutoff_    (para.cutoff    ()),
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

    std::string name() const noexcept override
    {
        return this->name_ + "BasePair "
               "(" + base_pair_.first + "-" + base_pair_.second + ")";
    }

  private:
    std::vector<indices_type> donor_indices_vec_;
    std::vector<indices_type> acceptor_indices_vec_;
    std::pair<std::string, std::string> base_pair_;
    index_pairs_type          ignore_list_;
    bool                      use_periodic_;
    std::string               ffgen_id_;

    double                        cutoff_;
    double                        alpha_BP_;
    double                        K_BP_;
    std::map<std::string, double> epsilon_BP_;
    std::map<std::string, double> r0_;
    std::map<std::string, double> theta0_1_;
    std::map<std::string, double> theta0_2_;
    std::map<std::string, double> phi0_;
    std::string                   name_;
};

#endif // OPEN_AICG2_PLUS_3SPN2_BASE_PAIR_FORCE_FIELD_GENERATOR_HPP
