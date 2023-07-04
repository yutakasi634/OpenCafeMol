#ifndef OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <OpenMM.h>
#include "ForceFieldGeneratorBase.hpp"
#include "src/util/Constants.hpp"

class ThreeSPN2ExcludedVolumeForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    ThreeSPN2ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& ignore_list, const bool use_periodic,
        const std::size_t ffgen_count = 0)
        : eps_(eps), cutoff_(cutoff), radiuses_(radiuses), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_str_(std::to_string(ffgen_count))
    {
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        std::string potential_formula = "energy;"
            "energy  = (epsilon*((sigma/r)^12 - 2*(sigma/r)^6) + epsilon) * step(sigma - r);"
            "epsilon = sqrt(epsilon1 * epsilon2);"
            "sigma   = 0.5*(  sigma1 +   sigma2);";

        const std::map<std::string, std::string> ff_params =
        {
            {"epsilon",  "TSPN2_EXV" + ffgen_id_str_ + "_epsilon"},
            {"sigma",    "TSPN2_EXV" + ffgen_id_str_ + "_sigma"},
        };

        for (auto itr=ff_params.begin(); itr!=ff_params.end(); ++itr)
        {
            potential_formula = std::regex_replace(
                potential_formula, std::regex(itr->first), itr->second);
        }

        auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        exv_ff->addPerParticleParameter(ff_params.at("epsilon"));
        exv_ff->addPerParticleParameter(ff_params.at("sigma"));

        double max_radius        = std::numeric_limits<double>::min();
        double second_max_radius = std::numeric_limits<double>::min();
        for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
        {
            const std::optional<double>& radius = radiuses_[idx];
            if(radius)
            {
                exv_ff->addParticle({eps_, radius.value()});

                if(max_radius < radius.value())
                {
                    second_max_radius = max_radius;
                    max_radius        = radius.value();
                }
            }
            else
            {
                exv_ff->addParticle({
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),});
            }
        }

        // if interaction_groups size is 0, no interaction group will be added,
        // so all the particle in the system will be considerd as participant
        for(const auto& group_pair : interaction_groups_)
        {
            exv_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance   = (max_radius + second_max_radius)*cutoff_;
        std::cerr << "        cutoff disntace is " << cutoff_distance << " nm" << std::endl;
        exv_ff->setCutoffDistance(cutoff_distance);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            exv_ff->addExclusion(pair.first, pair.second);
        }

        return exv_ff;
    }

    std::string name() const noexcept { return "3SPN2ExcludedVolume"; }

  private:
    double                              eps_;
    double                              cutoff_;
    std::vector<std::optional<double>>  radiuses_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
    bool                                use_periodic_;
    std::string                         ffgen_id_str_;
};


// ----------------------------------------------------------------------------
// 3SPN2 parameter set
//
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
//   - TABLE III (sigma)
//   - TABLE IV  (epsilon)

template<typename realT>
struct ThreeSPN2ExcludedVolumePotentialParameter
{
    using real_type = realT;
    real_type epsilon = real_type(1.0); // [kJ/mol]

    const std::map<std::string, real_type> sigma = { // [nm]
        {"P", real_type(4.5) * OpenMM::NmPerAngstrom},
        {"S", real_type(6.2) * OpenMM::NmPerAngstrom},
        {"A", real_type(5.4) * OpenMM::NmPerAngstrom},
        {"T", real_type(7.1) * OpenMM::NmPerAngstrom},
        {"G", real_type(4.9) * OpenMM::NmPerAngstrom},
        {"C", real_type(6.4) * OpenMM::NmPerAngstrom},
    };
};

#endif // OPEN_AICG2_PLUS_3SPN2_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
