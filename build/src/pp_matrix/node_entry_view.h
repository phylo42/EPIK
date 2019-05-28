#ifndef RAPPAS_CPP_NODE_ENTRY_VIEW_H
#define RAPPAS_CPP_NODE_ENTRY_VIEW_H

#include "row.h"
#include <boost/version.hpp>

/// static_vector has been accepted in boost only in v1.56. We check the version of boost library,
/// and use std::vector if the library is too old.
/// For our purposes, boost::static_vector is surely preferable due to its better performance.
#if ((BOOST_VERSION / 100 % 1000) < 56)
#define USE_NONSTATIC_VECTOR 1
#endif

#ifdef USE_NONSTATIC_VECTOR
#include <vector>
#pragma message("Boost is too old and does not provide static_vector. Using std::vector instead.")

namespace rappas
{
    namespace impl
    {
        template <class... Args>
        using stack_type = std::vector<Args...>;
    }
}
#else

#include <boost/container/static_vector.hpp>

namespace rappas
{
    namespace impl
    {
        template <typename T>
        using stack_type = boost::container::static_vector<T, core::seq_traits::max_kmer_length>;
    }
}
#endif


class node_entry;

/// \brief A forward access const iterator for phylo_kmer pairs [kmer value, score]. Iterates over
/// a fixed node_entry_view of size K.
class phylo_kmer_iterator
{
public:
    /// \brief Temporary data storage for mmers, m <= k. Used only in the branch-and-bound algorithm
    /// to generate phylo-kmers
    struct phylo_mmer
    {
        core::phylo_kmer mmer;
        /// a position of the last letter
        core::phylo_kmer::pos_type last_position;
        size_t last_index;
        size_t next_index;
    };

    /// Member types

    using iterator_category = std::forward_iterator_tag;
    using reference = const core::phylo_kmer&;
    using pointer = const core::phylo_kmer*;

    /// A stack type. Boost::static_vector of fixed size (core::seq_traits::max_kmer_length) if
    /// presented in boost, std::vector instead.
    using stack_type = rappas::impl::stack_type<phylo_mmer>;


    phylo_kmer_iterator(const node_entry* entry, size_t kmer_size,
                        core::phylo_kmer::pos_type start_pos, stack_type&& stack) noexcept;
    phylo_kmer_iterator(const phylo_kmer_iterator&) = delete;
    phylo_kmer_iterator(phylo_kmer_iterator&&) = default;
    phylo_kmer_iterator& operator=(const phylo_kmer_iterator& rhs) = delete;
    phylo_kmer_iterator& operator=(phylo_kmer_iterator&&) = delete;
    ~phylo_kmer_iterator() noexcept = default;

    bool operator==(const phylo_kmer_iterator& rhs) const noexcept;
    bool operator!=(const phylo_kmer_iterator& rhs) const noexcept;
    phylo_kmer_iterator& operator++();

    reference operator*();
    pointer operator->();

private:
    phylo_mmer next_phylokmer();

    const node_entry* _entry;
    size_t _kmer_size;
    const core::phylo_kmer::pos_type _start_pos;
    const core::phylo_kmer::score_type _threshold;
    stack_type _stack;
    phylo_mmer _current;
};

/// \brief A lightweight view of node_entry. Implements a "window" of size K over a node_entry.
class node_entry_view final
{
public:
    using const_iterator = phylo_kmer_iterator;
    using const_reference = const_iterator::reference;

    node_entry_view(const node_entry* entry, core::phylo_kmer::pos_type start, core::phylo_kmer::pos_type end) noexcept;
    node_entry_view(const node_entry_view& other) noexcept;
    node_entry_view(node_entry_view&&) = delete;
    node_entry_view& operator=(const node_entry_view&) = delete;
    node_entry_view& operator=(node_entry_view&& other) noexcept;
    ~node_entry_view() noexcept = default;

    const_iterator begin() const;
    const_iterator end() const;

    const node_entry* get_entry() const;
    core::phylo_kmer::pos_type get_start_pos() const;
    core::phylo_kmer::pos_type get_end_pos() const;

private:
    const node_entry* _entry;
    core::phylo_kmer::pos_type _start;
    core::phylo_kmer::pos_type _end;
};

bool operator==(const node_entry_view& a, const node_entry_view& b) noexcept;
bool operator!=(const node_entry_view& a, const node_entry_view& b) noexcept;


#endif //RAPPAS_CPP_NODE_ENTRY_VIEW_H
