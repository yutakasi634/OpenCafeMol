#ifndef OPEN_AICG2_PLUS_TOPOLOGY_HPP
#define OPEN_AICG2_PLUS_TOPOLOGY_HPP

#include <array>
#include <string>
#include <utility>
#include <vector>

class Topology
{
  public:
    using edge_type        = std::pair<std::size_t, std::size_t>;
    using index_pairs_type = std::vector<edge_type>;
    using molecule_type    = std::vector<std::size_t>;
    using node_type        = std::vector<std::pair<std::size_t, std::string>>;

  public:

    Topology(const std::size_t system_size)
        : nodes_({}), molecule_updated_(false)
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

    void make_molecule(const std::string& edge_type);

    index_pairs_type ignore_list_within_molecule() const;

    // This function return the vector of index pair in which first < second.
    index_pairs_type ignore_list_within_edge(
            const std::size_t dist, const std::string& edge_type) const;

    std::size_t                   size()  const noexcept { return nodes_.size(); }
    std::vector<node_type> const& nodes() const noexcept { return nodes_; }

  private:

    void list_adjacent_within(
        const std::size_t node_idx, const std::size_t dist,
        const std::string& edge_type, std::vector<std::size_t>& retval) const;

  private:
    std::vector<node_type>     nodes_;
    std::vector<molecule_type> molecules_;
    bool                       molecule_updated_;
};

#endif // OPEN_AICG2_PLUS_TOPOLOGY_HPP
