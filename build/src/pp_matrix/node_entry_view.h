#ifndef RAPPAS_CPP_NODE_ENTRY_VIEW_H
#define RAPPAS_CPP_NODE_ENTRY_VIEW_H

#include "row.h"

namespace rappas
{
    namespace impl
    {
        template<class T>
        using vector_type = std::vector<T>;
    }
}

class node_entry;

namespace rappas
{
    namespace impl
    {
        /// \brief Temporary data storage for mmers, m <= k. Used only in the branch-and-bound algorithm
        /// to generate phylo-kmers
        struct phylo_mmer
        {
            xpas::unpositioned_phylo_kmer mmer;
            /// a position of the last letter
            xpas::phylo_kmer::pos_type last_position;
            size_t last_index;
            size_t next_index;
        };

        /// \brief Divide-and-conquer phylo-kmer iterator.
        class dac_kmer_iterator
        {
        public:
            /// Member types
            using iterator_category = std::forward_iterator_tag;
            using reference = const xpas::unpositioned_phylo_kmer&;
            using pointer = const xpas::unpositioned_phylo_kmer*;

            dac_kmer_iterator(const node_entry* entry, size_t kmer_size, xpas::phylo_kmer::score_type threshold,
                              xpas::phylo_kmer::pos_type start_pos) noexcept;
            dac_kmer_iterator(const dac_kmer_iterator&) = delete;
            dac_kmer_iterator(dac_kmer_iterator&&) = default;
            dac_kmer_iterator& operator=(const dac_kmer_iterator&) = delete;
            dac_kmer_iterator& operator=(dac_kmer_iterator&& rhs) noexcept;
            ~dac_kmer_iterator() noexcept = default;

            bool operator==(const dac_kmer_iterator& rhs) const noexcept;
            bool operator!=(const dac_kmer_iterator& rhs) const noexcept;
            dac_kmer_iterator& operator++();

            reference operator*() const noexcept;
            pointer operator->() const noexcept;

        private:
            xpas::unpositioned_phylo_kmer _next_phylokmer();
            void _select_right_halfmers_bound();

            const node_entry* _entry;
            size_t _kmer_size;
            size_t _left_part_size;
            xpas::phylo_kmer::pos_type _start_pos;
            xpas::phylo_kmer::score_type _threshold;
            xpas::unpositioned_phylo_kmer _current;

            vector_type<xpas::unpositioned_phylo_kmer> _left_halfmers;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _left_halfmer_it;

            vector_type<xpas::unpositioned_phylo_kmer> _right_halfmers;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _right_halfmer_it;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _last_right_halfmer_it;
        };
    }
}


/// \brief A lightweight view of node_entry. Implements a "window" of size K over a node_entry.
class node_entry_view final
{
public:
    using const_iterator = rappas::impl::dac_kmer_iterator;
    using const_reference = const_iterator::reference;

    node_entry_view(const node_entry* entry, xpas::phylo_kmer::score_type threshold,
                    xpas::phylo_kmer::pos_type start, xpas::phylo_kmer::pos_type end) noexcept;
    node_entry_view(const node_entry_view& other) noexcept;
    node_entry_view(node_entry_view&&) = delete;
    node_entry_view& operator=(const node_entry_view&) = delete;
    node_entry_view& operator=(node_entry_view&& other) noexcept;
    ~node_entry_view() noexcept = default;

    const_iterator begin() const;
    const_iterator end() const noexcept;

    const node_entry* get_entry() const noexcept;
    xpas::phylo_kmer::pos_type get_start_pos() const noexcept;
    xpas::phylo_kmer::pos_type get_end_pos() const noexcept;
    xpas::phylo_kmer::score_type get_threshold() const noexcept;

private:
    const node_entry* _entry;
    xpas::phylo_kmer::score_type _threshold;
    xpas::phylo_kmer::pos_type _start;
    xpas::phylo_kmer::pos_type _end;
};

bool operator==(const node_entry_view& a, const node_entry_view& b) noexcept;
bool operator!=(const node_entry_view& a, const node_entry_view& b) noexcept;


#endif //RAPPAS_CPP_NODE_ENTRY_VIEW_H
