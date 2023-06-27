#ifndef OPEN_AICG2_PLUS_TOPOLOGY_HPP
#define OPEN_AICG2_PLUS_TOPOLOGY_HPP

#include <algorithm>
#include <limits>

class Topology
{
  public:
    using edge_type        = std::pair<std::size_t, std::size_t>;
    using index_pairs_type = std::vector<edge_type>;
    using molecule_type    = std::vector<std::size_t>;
    using node_type        = std::vector<std::pair<std::size_t, std::string>>;

  public:
    Topology(const std::size_t system_size) : nodes_({}), molecule_updated_(false)
    {
        nodes_.resize(system_size);
    }

    void add_edges(const std::vector<edge_type>& edges,
                   const std::string& edge_type)
    {
        for(const auto& edge : edges)
        {
            const std::size_t first  = edge.first;
            const std::size_t second = edge.second;

            nodes_[first] .push_back(std::make_pair(second, edge_type));
            nodes_[second].push_back(std::make_pair(first,  edge_type));
        }

        molecule_updated_ = false;
    }

    template<long unsigned int elemNum>
    void add_edges(const std::vector<std::array<std::size_t, elemNum>>& index_set_vec,
                   const std::string& edge_type)
    {
        for(const auto& indices : index_set_vec)
        {
            for(std::size_t idx_i=0; idx_i<indices.size()-1; ++idx_i)
            {
                for(std::size_t idx_j=idx_i+1; idx_j<indices.size(); ++idx_j)
                {
                    const std::size_t first  = indices[idx_i];
                    const std::size_t second = indices[idx_j];

                    nodes_[first] .push_back(std::make_pair(second, edge_type));
                    nodes_[second].push_back(std::make_pair(first,  edge_type));
                }
            }
        }
    }

    void make_molecule(const std::string& edge_type)
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

    index_pairs_type ignore_list_within_molecule() const
    {
        if(!molecule_updated_)
        {
            throw std::runtime_error(
                "[error] function ignore_list_within_molecule was called before "
                "updating molecule based on latest bond information by make_molecule.");
        }

        std::cerr << "        generating exclusion list from molecule information"
                  << std::endl;

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

    // This function return the vector of index pair in which first < second.
    index_pairs_type ignore_list_within_edge(
            const std::size_t dist, const std::string& edge_type) const
    {

        std::cerr << "        generating exclusion list from " << edge_type
                  << " information" << std::endl;

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

    std::size_t                   size()  const { return nodes_.size(); }
    const std::vector<node_type>& nodes() const { return nodes_; }

  private:
    void list_adjacent_within(
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

  private:
    std::vector<node_type>     nodes_;
    std::vector<molecule_type> molecules_;
    bool                       molecule_updated_;
};

#endif // OPEN_AICG2_PLUS_TOPOLOGY_HPP
