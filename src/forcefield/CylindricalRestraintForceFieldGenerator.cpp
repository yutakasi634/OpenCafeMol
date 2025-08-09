#include "CylindricalRestraintForceFieldGenerator.hpp"
#include "ForceFieldIDGenerator.hpp"

std::unique_ptr<OpenMM::Force> CylindricalRestraintForceFieldGenerator::generate() const
{
	std::string potential_formula = fmt::format(
		"{id}_k * (dr - {id}_v0)^2;"
		"dr = ((x - xh)^2 + (y - yh)^2 + (z - zh)^2)^(1/2);"
		"xh = ra_r * {id}_xa + {id}_x0;"
		"yh = ra_r * {id}_ya + {id}_y0;"
		"zh = ra_r * {id}_za + {id}_z0;"
		"ra_r = {id}_xa * x + {id}_ya * y + {id}_za * z;",
		fmt::arg("id", ffgen_id_));
	if (use_periodic_)
		potential_formula = fmt::format(
			"{id}_k * (dr - {id}_v0)^2;"
			"dr = periodicdistance(x, y, z, xh, yh, zh);"
			"xh = ra_r * {id}_xa + {id}_x0;"
			"yh = ra_r * {id}_ya + {id}_y0;"
			"zh = ra_r * {id}_za + {id}_z0;"
			"ra_r = {id}_xa * x + {id}_ya * y + {id}_za * z;",
		fmt::arg("id", ffgen_id_));

	auto external_ff = std::make_unique<OpenMM::CustomExternalForce>(potential_formula);
	external_ff->addPerParticleParameter(fmt::format("{}_xa", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_ya", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_za", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_x0", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_y0", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_z0", ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_k" , ffgen_id_));
	external_ff->addPerParticleParameter(fmt::format("{}_v0", ffgen_id_));

	for (std::size_t idx = 0; idx < indices_.size(); ++idx)
	{
		const auto axis = axes_[idx];
		const auto r0   = shifts_[idx];
		const double k  = ks_[idx];
		const double v0 = v0s_[idx];
		external_ff->addParticle(indices_[idx],
			{axis[0], axis[1], axis[2], r0[0], r0[1], r0[2], k, v0});
	}
	return external_ff;
}
