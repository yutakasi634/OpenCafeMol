#ifndef OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "../util/Constants.hpp"

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
          }, use_periodic_(use_periodic)
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

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        auto spline_func = std::make_unique<OpenMM::Continuous1DFunction>(
                std::vector<double>(spline_table_.begin(), spline_table_.end()),
                min_theta_, max_theta_);
        const double min_theta_y = spline_table_[0];
        const double max_theta_y = spline_table_[9];

        // calculate minimum energy of this potential
        double min_energy = min_theta_y;
        double th         = min_theta_;
        while(th < max_theta_)
        {
            const double energy = this->spline_interpolate(th);
            min_energy = std::min(min_energy, energy);
            th += 1e-4;
        }

        const std::string potential_formula =
            "k * ("
                "spline(theta) - min_energy"
                "+ (step(min_theta-theta) * (30 * (min_theta-theta) + min_theta_y))"
                "+ (step(theta-max_theta) * (30 * (theta-max_theta) + max_theta_y))"
            ");"
            "theta = angle(p1,p2,p3)";
        auto angle_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(3, potential_formula);
        angle_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
        angle_ff->addTabulatedFunction("spline", spline_func.release());
        angle_ff->addPerBondParameter("k");
        angle_ff->addGlobalParameter("min_energy",  min_energy);
        angle_ff->addGlobalParameter("min_theta",   min_theta_);
        angle_ff->addGlobalParameter("max_theta",   max_theta_);
        angle_ff->addGlobalParameter("min_theta_y", min_theta_y);
        angle_ff->addGlobalParameter("max_theta_y", max_theta_y);


        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const std::array<std::size_t, 3>& triplet = indices_vec_[idx];
            angle_ff->addBond({triplet.begin(), triplet.end()}, {ks_[idx]});
        }

        return angle_ff;
    }

    void add_exclusion(std::vector<std::pair<std::size_t, std::size_t>>& exclusion_pairs) const noexcept
    {
        for(const auto& indices : indices_vec_)
        {
            // TODO
            // dupulication in exclusion list make error.
            // all exclusion shoul be specified in HarmonicBond and GoContact

            exclusion_pairs.push_back(std::make_pair(indices[0], indices[2]));
        }
    }

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

  private:
    double spline_interpolate(const double th) const noexcept
    {
        constexpr double one_over_six = double(1.0) / double(6.0);

        std::size_t n = std::min<std::size_t>(9, std::floor((th - min_theta_) * rdtheta_));

        if      (thetas_[n+1] < th) { n++; }
        else if (th < thetas_[n])   { n--; }
        assert(n <= 9);
        assert(thetas_[n] <= th && th <= thetas_[n+1]);

        const double a = (thetas_[n+1] - th) * rdtheta_;
        const double b = (th - thetas_[n  ]) * rdtheta_;

        const double e1 = a * spline_table_[n] + b * spline_table_[n+1];
        const double e2 =
        ((a * a * a - a) * spline_second_deriv_table_[n] +
         (b * b * b - b) * spline_second_deriv_table_[n+1]) *
        dtheta_ * dtheta_ * one_over_six;

        return e1 + e2;
   }

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
    const bool                use_periodic_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
