#pragma once

#include <string>
#include <vector>
#include <unordered_map>

class phyml_result_parser;
class phylo_tree;

/// \brief A posterior probabilities matrix class.
/// \details A matrix class for storing posterior probabilities, given by the ancestral reconstruction
/// algorithm. So, this is a matrix of size [#branch_nodes x #sites x #variants], where:
/// - #branch_nodes is the number of non-leaf nodes of input tree
/// - #sites is the size of input alignment,
/// - #variants is the alphabet size.
class proba_matrix
{
    friend phyml_result_parser;
public:
    /// a row of size #variants, e.g. [0.1 0.3 0.7. 0.1]
    using pos_probs_t = std::vector<float>;
    /// a collection of rows for all the alignment positions
    using branch_entry_t = std::vector<pos_probs_t>;

    proba_matrix();
    proba_matrix(const proba_matrix&) = delete;
    proba_matrix(proba_matrix&& other);
    ~proba_matrix() = default;

    proba_matrix& operator=(const proba_matrix&) = delete;

    size_t num_branches() const;
    size_t num_sites() const;
    size_t num_variants() const;

    const branch_entry_t at(int branch_id) const;

private:
    void _add_branch_entry(int node_id, const branch_entry_t& branch_entry);

private:
    std::unordered_map<int, branch_entry_t> _data;
};
