#ifndef RAPPAS_CPP_BRANCH_ENTRY_VIEW_H
#define RAPPAS_CPP_BRANCH_ENTRY_VIEW_H

#include "row.h"
#include "plf_stack.h"

class branch_entry;


/// \brief A forward access const iterator for phylo_kmer pairs [kmer value, score]. Iterates over
/// a fixed branch_entry_view of size K.
class phylo_kmer_iterator
{
public:
    /// Temporary data storage for mmers, m <= k
    struct phylo_mmer
    {
        phylo_kmer mmer;
        size_t last_position;
        size_t last_index;
        bool visited;
    };

    using iterator_category = std::forward_iterator_tag;
    using reference = const phylo_kmer&;
    using pointer = const phylo_kmer*;
    using stack_type = plf::stack<phylo_mmer>;

    phylo_kmer_iterator(const branch_entry* entry, size_t kmer_size,
                        size_t start_pos, stack_type stack) noexcept;
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

    const branch_entry* _entry;
    size_t _kmer_size;
    size_t _start_pos;
    stack_type _stack;
    const score_t _threshold;
    phylo_mmer _current;
};

/// \brief A lightweight view of branch_entry. Implements a "window" of size K over a branch_entry.
class branch_entry_view final
{
public:
    using const_iterator = phylo_kmer_iterator;
    using const_reference = const_iterator::reference;

    branch_entry_view(const branch_entry* entry, size_t start, size_t end) noexcept;
    branch_entry_view(const branch_entry_view& other) noexcept;
    branch_entry_view(branch_entry_view&&) = delete;
    branch_entry_view& operator=(const branch_entry_view&) = delete;
    branch_entry_view& operator=(branch_entry_view&& other) noexcept;
    ~branch_entry_view() noexcept = default;

    const_iterator begin() const;
    const_iterator end() const;

    const branch_entry* get_entry() const;
    size_t get_start_pos() const;
    size_t get_end_pos() const;

private:
    const branch_entry* _entry;
    size_t _start;
    size_t _end;
};

bool operator==(const branch_entry_view& a, const branch_entry_view& b) noexcept;
bool operator!=(const branch_entry_view& a, const branch_entry_view& b) noexcept;


#endif //RAPPAS_CPP_BRANCH_ENTRY_VIEW_H
