#ifndef RAPPAS_CPP_PROBA_MATRIX_H
#define RAPPAS_CPP_PROBA_MATRIX_H

#include <unordered_map>
#include "node_entry.h"

namespace rappas
{

    /// \brief A posterior probabilities matrix class.
    /// \details A matrix class for storing posterior probabilities, given by the ancestral reconstruction
    /// algorithm. So, this is a matrix of size [#branch_nodes x #sites x #variants], where:
    /// - #branch_nodes is the number of non-leaf nodes of input tree
    /// - #sites is the size of input alignment,
    /// - #variants is the alphabet size.
    class proba_matrix final
    {
    public:
        static const branch_type NOT_A_LABEL = std::numeric_limits<branch_type>::max();

        /// a map for a fast access to a submatrix by branch node label
        using storage = std::unordered_map<std::string, node_entry>;
        using iterator = typename storage::iterator;
        using const_iterator = typename storage::const_iterator;
        using mapped_type = storage::mapped_type;

        proba_matrix() = default;
        proba_matrix(const proba_matrix&) = delete;
        proba_matrix(proba_matrix&& other) = default;
        proba_matrix& operator=(const proba_matrix&) = delete;
        proba_matrix& operator=(proba_matrix&&) = delete;
        ~proba_matrix() = default;

        /// capacity
        size_t num_branches() const;
        size_t num_sites() const;

        // Lookup
        mapped_type& operator[](const std::string& ar_label);
        const mapped_type& at(const std::string& ar_label) const;
        iterator find(const std::string& ar_label);
        const_iterator find(const std::string& ar_label) const;

        /// iterators
        iterator begin();
        iterator end();
        const_iterator begin() const;
        const_iterator end() const;

    private:
        storage _data;
    };
}

#endif