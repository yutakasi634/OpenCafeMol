#include "Topology.hpp"

#include "util/Logger.hpp"

#include <algorithm>

void Topology::make_molecule(const std::string& edge_type)
{
    molecules_.clear();
    std::vector<bool> node_anotated(nodes_.size(), false);

    const auto& begin_itr = node_anotated.begin();
    const auto& end_itr   = node_anotated.end();

    std::vector<std::size_t> search_stack;
    while(true)
    {
        const auto& unanotated_idx_itr = std::find(begin_itr, end_itr, false);

        if(unanotated_idx_itr == end_itr) { break; }

        // make new molecule
        const std::size_t unanotated_idx = std::distance(begin_itr, unanotated_idx_itr);
        search_stack.push_back(unanotated_idx);
        molecules_.push_back(molecule_type());
        molecule_type& molecule = molecules_.back();

        // construct member of molecule
        while(!search_stack.empty())
        {
            const std::size_t subject_idx = search_stack.back();
            molecule.push_back(subject_idx);
            node_anotated[subject_idx] = true;
            search_stack.pop_back();

            for(const auto& edge : nodes_[subject_idx])
            {
                const std::size_t adjacent_idx = edge.first;
                if(edge.second == edge_type && !node_anotated[adjacent_idx])
                {
                    search_stack.push_back(adjacent_idx);
                }
            }
        }
    }

    molecule_updated_ = true;
}


Topology::index_pairs_type Topology::ignore_list_within_molecule() const
{
    log_assert(molecule_updated_,
        "function ignore_list_within_molecule was called before "
        "updating molecule based on latest bond information by make_molecule.");

    log_info("        generating exclusion list from molecule information");

    std::vector<std::pair<std::size_t, std::size_t>> ignore_list;

    for(const auto& molecule : molecules_)
    {
        for(std::size_t idx_1st=0; idx_1st<molecule.size()-1; ++idx_1st)
        {
            for(std::size_t idx_2nd=idx_1st+1; idx_2nd<molecule.size(); ++idx_2nd)
            {
                ignore_list.push_back({molecule[idx_1st], molecule[idx_2nd]});
            }
        }
    }

    return ignore_list;
}

Topology::index_pairs_type Topology::ignore_list_within_edge(
        const std::size_t dist, const std::string& edge_type) const
{
    log_info("        generating exclusion list from {} information", edge_type);

    index_pairs_type ignore_list;
    for(std::size_t idx=0; idx<nodes_.size(); ++idx)
    {
        std::vector<std::size_t> within_bond_list;
        this->list_adjacent_within(idx, dist, edge_type, within_bond_list);
        for(std::size_t pair_idx : within_bond_list)
        {
            if(pair_idx != idx)
            {
                ignore_list.push_back({idx, pair_idx});
            }
        }
    }

    for(auto& pair : ignore_list)
    {
        if(pair.first > pair.second)
        {
            const std::size_t first  = pair.first;
            const std::size_t second = pair.second;
            pair.first = second;
            pair.second = first;
        }
    }

    std::sort(ignore_list.begin(), ignore_list.end());
    const auto& result = std::unique(ignore_list.begin(), ignore_list.end());
    ignore_list.erase(result, ignore_list.end());

    return ignore_list;
}



void Topology::list_adjacent_within(
    const std::size_t node_idx, const std::size_t dist,
    const std::string& edge_type, std::vector<std::size_t>& retval) const
{
    retval.push_back(node_idx);
    if(dist == 0) { return; }

    for(auto edge : nodes_[node_idx])
    {
        if(edge.second == edge_type)
        {
            const std::size_t edge_idx = edge.first;
            const bool not_added =
                (std::find(retval.begin(), retval.end(), edge_idx) == retval.end());
            if(not_added)
            {
                this->list_adjacent_within(edge_idx, dist-1, edge_type, retval);
            }
        }
    }
    return;
}
