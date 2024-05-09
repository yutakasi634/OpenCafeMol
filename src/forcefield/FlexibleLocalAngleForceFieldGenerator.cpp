#include "FlexibleLocalAngleForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> FlexibleLocalAngleForceFieldGenerator::generate() const
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

    const std::string potential_formula = fmt::format("{id}_k * ("
            "spline(theta) - {id}_min_energy + "
            "(step({id}_min_theta-theta) * (30 * ({id}_min_theta-theta) + {id}_min_theta_y)) +"
            "(step(theta-{id}_max_theta) * (30 * (theta-{id}_max_theta) + {id}_max_theta_y))"
        "); theta = angle(p1, p2, p3)", fmt::arg("id", ffgen_id_));

    auto angle_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(3, potential_formula);
    angle_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
    angle_ff->addTabulatedFunction("spline", spline_func.release());
    angle_ff->addPerBondParameter(fmt::format("{}_k",           ffgen_id_));
    angle_ff->addGlobalParameter (fmt::format("{}_min_energy",  ffgen_id_), min_energy);
    angle_ff->addGlobalParameter (fmt::format("{}_min_theta",   ffgen_id_), min_theta_);
    angle_ff->addGlobalParameter (fmt::format("{}_max_theta",   ffgen_id_), max_theta_);
    angle_ff->addGlobalParameter (fmt::format("{}_min_theta_y", ffgen_id_), min_theta_y);
    angle_ff->addGlobalParameter (fmt::format("{}_max_theta_y", ffgen_id_), max_theta_y);


    for(std::size_t idx=0; idx<indices_vec_.size(); ++idx)
    {
        const std::array<std::size_t, 3>& triplet = indices_vec_[idx];
        angle_ff->addBond({triplet.begin(), triplet.end()}, {ks_[idx]});
    }

    return angle_ff;
}

double FlexibleLocalAngleForceFieldGenerator::spline_interpolate(const double th) const noexcept
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


