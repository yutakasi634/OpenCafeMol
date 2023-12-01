#ifndef OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP

#include <OpenMM.h>
#include <fmt/core.h>
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

class UniformWeeksChandlerAndersenForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using index_pairs_type       = std::vector<std::pair<std::size_t, std::size_t>>;
    using interaction_group_type = std::pair<std::set<int>, std::set<int>>;

  public:
    UniformWeeksChandlerAndersenForceFieldGenerator(
        const std::size_t system_size, const double eps, const double sigma,
        const std::vector<std::size_t> former_group_vec,
        const std::vector<std::size_t> latter_group_vec,
        const index_pairs_type& ignore_list,
        const bool use_periodic,
        const std::vector<std::pair<std::string, std::string>> ignore_group_pairs = {},
        const std::vector<std::optional<std::string>> group_vec = {})
        : system_size_(system_size), eps_(eps), sigma_(sigma), ignore_list_(ignore_list),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("UWCA{}", ffid.gen())),
          former_group_size_(former_group_vec.size()), latter_group_size_(latter_group_vec.size())
    {
        // make interaction group
        if(ignore_group_pairs.size() == 0)
        {
            std::set<int> former_participants;
            for(const std::size_t idx : former_group_vec)
            {
                former_participants.insert(idx);
            }
            std::set<int> latter_participants;
            for(const std::size_t idx : latter_group_vec)
            {
                latter_participants.insert(idx);
            }
            interaction_groups_.push_back({ former_participants, latter_participants });
        }
        else // group based ignoration specified case
        {
            std::set<std::string> related_group_names;
            for(const auto& name_group_pair : ignore_group_pairs)
            {
                related_group_names.insert(name_group_pair.first);
                related_group_names.insert(name_group_pair.second);
            }

            std::map<std::string, std::set<int>> former_related_group_map;
            std::map<std::string, std::set<int>> latter_related_group_map;
            for(const auto& name : related_group_names)
            {
                former_related_group_map.insert(std::make_pair(name, std::set<int>()));
                latter_related_group_map.insert(std::make_pair(name, std::set<int>()));
            }

            std::set<int> former_others;
            for(const std::size_t idx : former_group_vec)
            {
                if(group_vec[idx])
                {
                    const std::string group_name = group_vec[idx].value();
                    if(related_group_names.count(group_name) != 0)
                    {
                        former_related_group_map.at(group_name).insert(idx);
                    }
                    else
                    {
                        former_others.insert(idx);
                    }
                }
                else
                {
                    former_others.insert(idx);
                }
            }

            std::set<int> latter_others;
            for(const std::size_t idx : latter_group_vec)
            {
                if(group_vec[idx])
                {
                    const std::string group_name = group_vec[idx].value();
                    if(related_group_names.count(group_name) != 0)
                    {
                        latter_related_group_map.at(group_name).insert(idx);
                    }
                    else
                    {
                        latter_others.insert(idx);
                    }
                }
                else
                {
                    latter_others.insert(idx);
                }
            }

            if(!former_others.empty())
            {
                for(const auto& latter_name_set : latter_related_group_map)
                {
                    const std::set<int>& second_group = latter_name_set.second;
                    if(second_group.empty()) { continue; }
                    interaction_groups_.push_back({ former_others, second_group });
                }

                if(!latter_others.empty())
                {
                    interaction_groups_.push_back({ former_others, latter_others });
                }
            }

            for(const auto& former_name_set : former_related_group_map)
            {
                const std::set<int>& first_group = former_name_set.second;
                if(first_group.empty()) { continue; }
                if(!latter_others.empty())
                {
                    interaction_groups_.push_back({ first_group, latter_others });
                }

                const std::string& first_name = former_name_set.first;
                for(const auto& latter_name_set : latter_related_group_map)
                {
                    const std::string& second_name = latter_name_set.first;
                    if(!Utility::contains(ignore_group_pairs, { first_name, second_name }))
                    {
                        const std::set<int>& second_group = latter_name_set.second;
                        if(!second_group.empty())
                        {
                            interaction_groups_.push_back({ first_group, second_group });
                        }
                    }
                }
            }
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const noexcept override
    {
        const std::string potential_formula = fmt::format(
            "4 * {id}_eps * (sigma_r_12 - sigma_r_6 + 0.25);"
            "sigma_r_12 = sigma_r_6^2;"
            "sigma_r_6  = sigma_r^6;"
            "sigma_r    = {id}_sigma / r",
            fmt::arg("id", ffgen_id_));

        auto uwca_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

        uwca_ff->addGlobalParameter(fmt::format("{}_sigma", ffgen_id_), sigma_);
        uwca_ff->addGlobalParameter(fmt::format("{}_eps", ffgen_id_), eps_);

        for(std::size_t idx=0; idx<system_size_; ++idx)
        {
            uwca_ff->addParticle();
        }

        for(const auto& group_pair : interaction_groups_)
        {
            uwca_ff->addInteractionGroup(group_pair.first, group_pair.second);
        }

        // set pbc condition
        if(use_periodic_)
        {
            uwca_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
        }
        else
        {
            uwca_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
        }

        // set cutoff
        const double cutoff_distance = sigma_ * std::pow(2.0, 1.0/6.0);
        std::cerr << "    UniformWeeksChandlerAndersen : cutoff disntace is "
                  << cutoff_distance << " nm" << std::endl;
        uwca_ff->setCutoffDistance(cutoff_distance);

        // set excludion list
        for(const auto& pair : ignore_list_)
        {
            uwca_ff->addExclusion(pair.first, pair.second);
        }

        return uwca_ff;
    }

    std::size_t former_group_size() const noexcept { return former_group_size_; }
    std::size_t latter_group_size() const noexcept { return latter_group_size_; }

    std::string name() const noexcept { return "UniformWeeksChadlerAndersen"; }

  private:
    std::size_t      system_size_;
    double           eps_;
    double           sigma_;
    index_pairs_type ignore_list_;
    bool             use_periodic_;
    std::string      ffgen_id_;
    std::size_t      former_group_size_;
    std::size_t      latter_group_size_;

    std::vector<interaction_group_type> interaction_groups_;
};

#endif // OPEN_AICG2_PLUS_UNIFORM_WEEKS_CHANDLER_ANDERSEN_FORCE_FIELD_GENERATOR_HPP
