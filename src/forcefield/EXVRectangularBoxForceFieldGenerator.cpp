#include "src/util/Logger.hpp"
#include "EXVRectangularBoxForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <iostream>

EXVRectangularBoxForceFieldGenerator::EXVRectangularBoxForceFieldGenerator(
    const double eps, const std::array<double, 3>& box_lower, const std::array<double, 3>& box_upper,
    const std::vector<std::size_t>& indices, const std::vector<double>& radii,
    const bool use_periodic)
    : eps_(eps), box_lower_(box_lower), box_upper_(box_upper), indices_(indices), radii_(radii),
      use_periodic_(use_periodic),
      ffgen_id_(fmt::format("EXVRB{}", ffid.gen()))
{
    if(use_periodic_)
    {
        throw std::runtime_error(
            "[error] read_toml_rectangular_box_ff_generator :"
            " RectangularBox interaction does not work under "
            "periodic boundary condition. Use unlimited boundary.");
    }

    const std::string log_template =
        "{} coordinate of the box's upper edge must be larger than that of the box's lower edge.";
    log_assert(box_upper[0] > box_lower[0], log_template, "X");
    log_assert(box_upper[1] > box_lower[1], log_template, "Y");
    log_assert(box_upper[2] > box_lower[2], log_template, "Z");
}

std::unique_ptr<OpenMM::Force> EXVRectangularBoxForceFieldGenerator::generate() const
{
    const std::string potential_formula = fmt::format(
        "{id}_epsilon * ("
        "    ({id}_radius/(x-{id}_x_lower))^12 + ({id}_radius/({id}_x_upper-x))^12 +"
        "    ({id}_radius/(y-{id}_y_lower))^12 + ({id}_radius/({id}_y_upper-y))^12 +"
        "    ({id}_radius/(z-{id}_z_lower))^12 + ({id}_radius/({id}_z_upper-z))^12)",
        fmt::arg("id", ffgen_id_));

    auto box_ff = std::make_unique<OpenMM::CustomExternalForce>(potential_formula);

    box_ff->addPerParticleParameter(fmt::format("{}_radius", ffgen_id_));
    box_ff->addGlobalParameter(fmt::format("{}_epsilon", ffgen_id_), eps_);
    box_ff->addGlobalParameter(fmt::format("{}_x_lower", ffgen_id_), box_lower_[0]);
    box_ff->addGlobalParameter(fmt::format("{}_y_lower", ffgen_id_), box_lower_[1]);
    box_ff->addGlobalParameter(fmt::format("{}_z_lower", ffgen_id_), box_lower_[2]);
    box_ff->addGlobalParameter(fmt::format("{}_x_upper", ffgen_id_), box_upper_[0]);
    box_ff->addGlobalParameter(fmt::format("{}_y_upper", ffgen_id_), box_upper_[1]);
    box_ff->addGlobalParameter(fmt::format("{}_z_upper", ffgen_id_), box_upper_[2]);

    for(std::size_t param_idx=0; param_idx<indices_.size(); ++param_idx)
    {
        const double radius = radii_[param_idx];
        box_ff->addParticle(indices_[param_idx], {radius});
    }

    return box_ff;
}
