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

namespace rappas
{
    namespace impl
    {
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

        /// \brief Divide-and-conquer phylo-kmer iterator.
        class dac_kmer_iterator
        {
        public:
            /// Member types
            using iterator_category = std::forward_iterator_tag;
            using reference = const core::phylo_kmer&;
            using pointer = const core::phylo_kmer*;

            dac_kmer_iterator(const node_entry* entry, size_t kmer_size, core::phylo_kmer::score_type threshold,
                              core::phylo_kmer::pos_type start_pos) noexcept;
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
            core::phylo_kmer _next_phylokmer();
            void _select_right_halfmers_bound();

            const node_entry* _entry;
            size_t _kmer_size;
            size_t _left_part_size;
            core::phylo_kmer::pos_type _start_pos;
            core::phylo_kmer::score_type _threshold;
            core::phylo_kmer _current;

            std::vector<core::phylo_kmer> _left_halfmers;
            std::vector<core::phylo_kmer>::iterator _left_halfmer_it;

            std::vector<core::phylo_kmer> _right_halfmers;
            std::vector<core::phylo_kmer>::iterator _right_halfmer_it;
            std::vector<core::phylo_kmer>::iterator _last_right_halfmer_it;
        };
    }
}


/// \brief A lightweight view of node_entry. Implements a "window" of size K over a node_entry.
class node_entry_view final
{
public:
    using const_iterator = rappas::impl::dac_kmer_iterator;
    using const_reference = const_iterator::reference;

    node_entry_view(const node_entry* entry, core::phylo_kmer::score_type threshold,
        core::phylo_kmer::pos_type start, core::phylo_kmer::pos_type end) noexcept;
    node_entry_view(const node_entry_view& other) noexcept;
    node_entry_view(node_entry_view&&) = delete;
    node_entry_view& operator=(const node_entry_view&) = delete;
    node_entry_view& operator=(node_entry_view&& other) noexcept;
    ~node_entry_view() noexcept = default;

    const_iterator begin() const;
    const_iterator end() const noexcept;

    const node_entry* get_entry() const noexcept;
    core::phylo_kmer::pos_type get_start_pos() const noexcept;
    core::phylo_kmer::pos_type get_end_pos() const noexcept;
    core::phylo_kmer::score_type get_threshold() const noexcept;

private:
    const node_entry* _entry;
    core::phylo_kmer::score_type _threshold;
    core::phylo_kmer::pos_type _start;
    core::phylo_kmer::pos_type _end;
};

bool operator==(const node_entry_view& a, const node_entry_view& b) noexcept;
bool operator!=(const node_entry_view& a, const node_entry_view& b) noexcept;


#endif //RAPPAS_CPP_NODE_ENTRY_VIEW_H
