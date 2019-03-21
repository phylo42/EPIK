#ifndef RAPPAS_CPP_PROBA_MATRIX_H
#define RAPPAS_CPP_PROBA_MATRIX_H


#include <vector>
#include <unordered_map>


/// \brief A posterior probabilities matrix class.
/// \details A matrix class for storing posterior probabilities, given by the ancestral reconstruction
/// algorithm. So, this is a matrix of size [#branch_nodes x #sites x #variants], where:
/// - #branch_nodes is the number of non-leaf nodes of input tree
/// - #sites is the size of input alignment,
/// - #variants is the alphabet size.
class proba_matrix
{
public:
    using key_t = uint8_t;
    static const key_t NOT_A_LABEL = std::numeric_limits<key_t>::max();

    /// a row of size #variants, e.g. [0.1 0.3 0.7. 0.1]
    using row_probs_t = std::vector<float>;
    /// a row of size #variants that stores the ordering of row_probs_t, e.g. [3, 0, 2, 1] (T, A, G, c)
    using row_pos_t = std::vector<unsigned char>;
    /// a row of matrix
    using row_t = std::pair<row_probs_t, row_pos_t>;
    /// a collection of rows for all the alignment positions
    using branch_entry_t = std::vector<row_t>;
    /// a map for a fast access to a submatrix by branch node label
    using storage = std::unordered_map<key_t, branch_entry_t>;
    using iterator = typename storage::iterator;
    using const_iterator = typename storage::const_iterator;

    proba_matrix() = default;
    proba_matrix(const proba_matrix&) = delete;
    proba_matrix(proba_matrix&& other) = default;
    ~proba_matrix() = default;

    proba_matrix& operator=(const proba_matrix&) = delete;
    proba_matrix&& operator=(proba_matrix&&) = delete;

    /// capacity
    size_t num_branches() const;
    size_t num_sites() const;
    size_t num_variants() const;

    // Lookup
    branch_entry_t& operator[](key_t branch_id);
    const branch_entry_t at(key_t branch_id) const;
    iterator find(const key_t& key);
    const_iterator find(const key_t& key) const;

    /// iterators
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    void sort();

private:
    storage _data;
};

#endif