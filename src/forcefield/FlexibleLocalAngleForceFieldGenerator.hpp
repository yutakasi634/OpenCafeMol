#ifndef OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP

#include "src/util/Constants.hpp"
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <sstream>
#include <string>

class FlexibleLocalAngleForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::array<std::size_t, 3>;

  public:
    FlexibleLocalAngleForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::array<double, 10>& spline_table, const std::string& aa_name,
        const bool use_periodic,
        const double min_theta = Constant::fla_spline_min_theta,
        const double max_theta = Constant::fla_spline_max_theta)
        : indices_vec_(indices_vec), ks_(ks), aa_name_(aa_name),
          min_theta_(min_theta), max_theta_(max_theta),
          dtheta_((max_theta - min_theta) / 9.0), rdtheta_(1.0 / dtheta_),
          spline_table_(spline_table),
          spline_second_deriv_table_(Constant::fla_spline_second_derivative_table.at(aa_name)),
          thetas_{
              {1.30900, 1.48353, 1.65806, 1.83260, 2.00713,
               2.18166, 2.35619, 2.53073, 2.70526, 2.87979}
          }, use_periodic_(use_periodic), ffgen_id_(fmt::format("FLA{}", ffid.gen()))
    {
        if(!(indices_vec.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] FlexibleLocalAngleForceFieldGenerator : "
                   "parameter number of"
                   "indices_vec (" << indices_vec.size() << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

    std::string name() const override
    {
        return "FlexibleLocalAngle (" + aa_name_ + ")";
    }

  private:

    double spline_interpolate(const double th) const noexcept;

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::string               aa_name_;
    double                    min_theta_ys_;
    double                    max_theta_ys_;
    double                    min_theta_;
    double                    max_theta_;
    double                    dtheta_;
    double                    rdtheta_;
    std::array<double, 10>    spline_table_;
    std::array<double, 10>    spline_second_deriv_table_;
    std::array<double, 10>    thetas_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
