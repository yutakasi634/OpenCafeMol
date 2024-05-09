#ifndef OPEN_AICG2_PLUS_INTERACTION_GROUP_HPP
#define OPEN_AICG2_PLUS_INTERACTION_GROUP_HPP

#include "src/util/Logger.hpp"
#include "src/util/Utility.hpp"

#include <algorithm>
#include <optional>
#include <map>
#include <set>
#include <utility>
#include <vector>

// without interaction group force calculation is faster than
// with that. So if all particle have the parameter for this
// calculation, force don't use interaction group.
template<typename T>
std::vector<std::pair<std::set<int>, std::set<int>>>
extract_interaction_group(
    const std::vector<std::optional<T>>&                    para,
    const std::vector<std::pair<std::string, std::string>>& ignore_group_pairs,
    const std::vector<std::optional<std::string>>&          group_vec)
{
    std::vector<std::pair<std::set<int>, std::set<int>>> interaction_groups;

    if(ignore_group_pairs.empty()) // no groups related to this FF
    {
        if(std::any_of(para.begin(), para.end(),
              [](const std::optional<T>& val) { return ! val.has_value(); }))
        {
            std::set<int> participants;
            for(std::size_t idx=0; idx<para.size(); ++idx)
            {
                if(para.at(idx).has_value())
                {
                    participants.insert(static_cast<int>(idx));
                }
            }
            interaction_groups.emplace_back(participants, participants);
        }
        else // we don't need to consider groups. do nothing.
        {
            log_info("        all particles are participants in this interaction");
        }
    }
    else // group based ignoration specified case
    {
        std::set<std::string> related_group_names;
        for(const auto& [group_i, group_j] : ignore_group_pairs)
        {
            related_group_names.insert(group_i);
            related_group_names.insert(group_j);
        }

        std::map<std::string, std::set<int>> related_group_map;
        for(const auto& name : related_group_names)
        {
            related_group_map.emplace(name, std::set<int>{});
        }

        std::set<int> others; // particles that does not belong to groups ignored
        for(std::size_t idx=0; idx<para.size(); ++idx)
        {
            if (para.at(idx).has_value()      && // the particle interacts with this
                group_vec.at(idx).has_value() && // the particle belongs to a group
                related_group_names.count(group_vec.at(idx).value()) != 0) // the group is in ignoring pair
            {
                related_group_map.at(group_vec.at(idx).value()).insert(idx);
            }
            else
            {
                others.insert(idx);
            }
        }

        const std::vector<std::pair<std::string, std::set<int>>>
            related_group_vec(related_group_map.cbegin(), related_group_map.cend());

        interaction_groups.emplace_back(others, others);

        for(std::size_t i=0; i<related_group_vec.size(); ++i)
        {
            const auto& [name_i, group_i] = related_group_vec.at(i);

            interaction_groups.emplace_back(group_i, others);

            for(std::size_t j=i; j<related_group_vec.size(); ++j)
            {
                const auto& [name_j, group_j] = related_group_vec.at(j);

                if(!Utility::contains(ignore_group_pairs, { name_i, name_j }))
                {
                    interaction_groups.emplace_back(group_i, group_j);
                }
            }
        }
    }
    return interaction_groups;
}

#endif// OPEN_AICG2_PLUS_INTERACTION_GROUP_HPP
