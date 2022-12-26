#ifndef OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP

#include <memory>
#include <limits>
#include <set>
#include <optional>
#include <OpenMM.h>

class ExcludedVolumeForceFieldGenerator final: public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    ExcludedVolumeForceFieldGenerator(const double eps, const double cutoff,
        const std::vector<std::optional<double>>& radiuses,
        const index_pairs_type& ignore_list,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : eps_(eps), cutoff_(cutoff), radiuses_(radiuses), ignore_list_(ignore_list)
    {
        // make interaction group
        if(ignore_group_pairs.size() == 0)
        {
            std::set<int> participants;
            for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
            {
                if(radiuses_[idx])
                {
                    participants.insert(idx);
                }
            }
            interaction_groups_.push_back({ participants, participants });
        }
        else // group based ignoration specified case
        {
            std::set<std::string> related_group_names;
            for(const auto& name_group_pair : ignore_group_pairs)
            {
                related_group_names.insert(name_group_pair.first);
                related_group_names.insert(name_group_pair.second);
            }

            std::map<std::string, std::set<int>> related_group_map;
            for(const auto& name : related_group_names)
            {
                related_group_map.insert(std::make_pair(name, std::set<int>()));
            }

            std::set<int> others;
            for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
            {
                if(radiuses_[idx])
                {
                    if(group_vec[idx])
                    {
                        const std::string group_name = group_vec[idx].value();
                        if(related_group_names.count(group_name) != 0)
                        {
                            related_group_map.at(group_name).insert(idx);
                        }
                        else
                        {
                            others.insert(idx);
                        }
                    }
                    else
                    {
                        others.insert(idx);
                    }
                }
            }

            std::vector<std::pair<std::string, std::set<int>>> related_group_vec;
            for(const auto& name_group_pair : related_group_map)
            {
                related_group_vec.push_back(name_group_pair);
            }

            interaction_groups_.push_back({ others, others });
            for(std::size_t idx_i=0 ; idx_i<related_group_vec.size(); ++idx_i)
            {
                const auto& name_group_pair_i = related_group_vec[idx_i];
                const std::string&   first_name  = name_group_pair_i.first;
                const std::set<int>& first_group = name_group_pair_i.second;
                interaction_groups_.push_back({ first_group, others });
                for(std::size_t idx_j=idx_i; idx_j<related_group_vec.size(); ++idx_j)
                {
                    const auto& name_group_pair_j = related_group_vec[idx_j];
                    const std::string&   second_name  = name_group_pair_j.first;
                    if(!Utility::contains(ignore_group_pairs, { first_name, second_name }))
                    {
                        const std::set<int>& second_group = name_group_pair_j.second;
                        interaction_groups_.push_back({ first_group, second_group });
                    }
                }
            }
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula =
            "epsilon*(sigma1+sigma2)^12*((1/r)^12-cutoff_correction)";
        auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        exv_ff->addPerParticleParameter("sigma");
        exv_ff->addGlobalParameter("epsilon", eps_);

        double max_radius        = std::numeric_limits<double>::min();
        double second_max_radius = std::numeric_limits<double>::min();
        std::set<int> participants;
        for(std::size_t idx=0; idx<radiuses_.size(); ++idx)
        {
            const std::optional<double>& radius = radiuses_[idx];
            if(radius)
            {
                exv_ff->addParticle({radius.value()});
                participants.insert(idx);

                if(max_radius < radius.value())
                {
                    second_max_radius = max_radius;
                    max_radius        = radius.value();
                }
            }
            else
            {
                exv_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
            }
        }

        for(const auto& group_pair : interaction_groups_)
        {
            exv_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set cutoff
        exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        const double cutoff_distance   = (max_radius + second_max_radius)*cutoff_;
        const double cutoff_correction = std::pow(1.0 / cutoff_distance, 12);
        exv_ff->setCutoffDistance(cutoff_distance);
        exv_ff->addGlobalParameter("cutoff_correction", cutoff_correction);

        // set exclusion list
        for(const auto& pair : ignore_list_)
        {
            exv_ff->addExclusion(pair.first, pair.second);
        }

        return exv_ff;
    }

    void add_exclusion(index_pairs_type exclusion_pairs) noexcept
    {

        for(auto& pair : exclusion_pairs)
        {
            if(pair.first > pair.second)
            {
                const std::size_t first  = pair.first;
                const std::size_t second = pair.second;
                pair.first = second;
                pair.second = first;
            }
        }

        for(const auto& pair : exclusion_pairs)
        {
            const auto result =
                std::find(ignore_list_.begin(), ignore_list_.end(), pair);
            if(result == ignore_list_.end())
            {
                ignore_list_.push_back(std::make_pair(pair.first, pair.second));
            }
        }
    }

  private:
    double                              eps_;
    double                              cutoff_;
    std::vector<std::optional<double>>  radiuses_;
    index_pairs_type                    ignore_list_;
    std::vector<interaction_group_type> interaction_groups_;
};

#endif // OPEN_AICG2_PLUS_EXCLUDED_VOLUME_FORCE_FIELD_GENERATOR_HPP
