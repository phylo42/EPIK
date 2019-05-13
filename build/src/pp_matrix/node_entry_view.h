#ifndef RAPPAS_CPP_NODE_ENTRY_VIEW_H
#define RAPPAS_CPP_NODE_ENTRY_VIEW_H

#include <boost/container/static_vector.hpp>
#include "row.h"

class node_entry;

/// \brief A forward access const iterator for phylo_kmer pairs [kmer value, score]. Iterates over
/// a fixed node_entry_view of size K.
class phylo_kmer_iterator
{
public:
    /// Temporary data storage for mmers, m <= k
    struct phylo_mmer
    {
        core::phylo_kmer mmer;
        core::phylo_kmer::pos_type last_position;
        size_t last_index;
        bool visited;
    };

    using iterator_category = std::forward_iterator_tag;
    using reference = const core::phylo_kmer&;
    using pointer = const core::phylo_kmer*;
    using stack_type = boost::container::static_vector<phylo_mmer, core::seq_traits::max_kmer_length>;

    phylo_kmer_iterator(const node_entry* entry, size_t kmer_size,
                        core::phylo_kmer::pos_type start_pos, stack_type stack) noexcept;
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
    void next_index();
    void next_position();
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
