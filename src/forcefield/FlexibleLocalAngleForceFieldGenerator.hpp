#ifndef OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "../util/Constants.hpp"

class FlexibleLocalAngleForceFieldGenerator
{
  public:
    using indices_type = std::vector<std::size_t>;

  public:
    FlexibleLocalAngleForceFieldGenerator(
        const std::vector<indices_type>& indices_vec, const std::vector<double>& ks,
        const std::array<double, 10>& spline_table, const std::string& aa_name,
        const double min_theta = Constant::fla_spline_min_theta,
        const double max_theta = Constant::fla_spline_max_theta)
        : indices_vec_(indices_vec), ks_(ks), aa_name_(aa_name),
          min_theta_(min_theta), max_theta_(max_theta), spline_table_(spline_table)
    {
        const std::size_t system_size = indices_vec.size();
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

    std::unique_ptr<OpenMM::CustomCompoundBondForce> generate() const noexcept
    {
        std::cerr << "    BondAngle     : FlexibleLocalAngle ("
                  << aa_name_ << ")" << std::endl;

        auto spline_func = std::make_unique<OpenMM::Continuous1DFunction>(
                std::vector<double>(spline_table_.begin(), spline_table_.end()),
                min_theta_, max_theta_);
        const double min_theta_y = spline_table_[0];
        const double max_theta_y = spline_table_[9];

        const std::string potential_formula =
            "k * ("
                "spline(theta)"
                "- (step(min_theta-theta) * 30 * (theta-min_theta) + min_theta_y)"
                "+ (step(theta-max_theta) * 30 * (theta-max_theta) + max_theta_y)"
            ");"
            "theta = angle(p1,p2,p3)";
        auto angle_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(3, potential_formula);
        angle_ff->addTabulatedFunction("spline", spline_func.release());
        angle_ff->addPerBondParameter("k");
        angle_ff->addPerBondParameter("min_theta");
        angle_ff->addPerBondParameter("min_theta_y");
        angle_ff->addPerBondParameter("max_theta");
        angle_ff->addPerBondParameter("max_theta_y");

        for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
        {
            const std::vector<std::size_t>& triplet = indices_vec_[idx];
            angle_ff->addBond({triplet.begin(), triplet.end()},
                              {ks_[idx], min_theta_, min_theta_y, max_theta_, max_theta_y});
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

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::string               aa_name_;
    double                    min_theta_ys_;
    double                    max_theta_ys_;
    double                    min_theta_;
    double                    max_theta_;
    std::array<double, 10>    spline_table_;
};

#endif // OPEN_AICG2_PLUS_FLEXIBLE_LOCAL_ANGLE_FORCE_FIELD_GENERATOR_HPP
