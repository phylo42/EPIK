#ifndef RAPPAS_CPP_BRANCH_ENTRY_H
#define RAPPAS_CPP_BRANCH_ENTRY_H

#include <vector>
#include <core/seq_traits.h>
#include <utils/meta.h>
#include "row.h"

class branch_entry_view;

/// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
class branch_entry final
{
    template <bool IsConst>
    class sliding_window_iterator
    {
        sliding_window_iterator(const branch_entry_view& view);
    };

public:
    using iterator = sliding_window_iterator<false>;
    using const_iterator = sliding_window_iterator<true>;

    explicit branch_entry(const seq_traits& traits = dna_seq_traits);
    branch_entry(std::vector<row>&& rows, const seq_traits& traits);
    branch_entry(const branch_entry&) = delete;
    branch_entry(branch_entry&&) = default;
    branch_entry& operator=(const branch_entry&) = delete;
    branch_entry& operator=(branch_entry&&) = default;
    ~branch_entry() noexcept = default;

    void push_back(row&& r);

    size_t alignment_size() const;
    size_t alphabet_size() const;

    const proba_pair& at(size_t position, size_t variant) const;

private:
    std::vector<row> _rows;
    seq_traits _traits;
};

/// \brief A lightweight view of branch_entry. Implements a "window" of size K over a branch_entry.
class branch_entry_view final
{
    /// \brief A forward access (non-)const iterator for phylo_kmer pairs [kmer value, score]. Iterates over
    /// a fixed branch_entry_view of size K.
    template <bool IsConst>
    class phylo_kmer_iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef typename choose<IsConst, const phylo_kmer&, phylo_kmer&>::type reference;
        typedef typename choose<IsConst, const phylo_kmer*, phylo_kmer*>::type pointer;

        phylo_kmer_iterator(const branch_entry_view& view, const phylo_kmer& kmer,
            size_t start_pos, size_t last_letter_idx)
            : _view(view)
            , _kmer(kmer)
            , _start_pos(start_pos)
            , _last_letter_idx(last_letter_idx)
        {}
        phylo_kmer_iterator(const phylo_kmer_iterator&) = default;
        phylo_kmer_iterator(phylo_kmer_iterator&&) = delete;
        phylo_kmer_iterator&& operator=(phylo_kmer_iterator&&) = delete;
        ~phylo_kmer_iterator() noexcept = default;

        phylo_kmer_iterator& operator=(const phylo_kmer_iterator& rhs) = delete;
        /*{
            if (*this != rhs)
            {
                _view = rhs._view;
                _kmer = rhs._kmer;
                _start_pos = rhs._start_pos;
                _last_letter_idx = rhs._last_letter_idx;
            }
            return *this;
        }*/

        bool operator==(const phylo_kmer_iterator& rhs) const noexcept
        {
            return _kmer == rhs._kmer;
        }

        bool operator!=(const phylo_kmer_iterator& rhs) const noexcept
        {
            return !(*this == rhs);
        }

        phylo_kmer_iterator& operator++()
        {
            auto letter = _view.at(_start_pos, _last_letter_idx).index;
            auto next = shift_append(_kmer.value, _view.kmer_size(), _view.alphabet_size(), letter);
            return phylo_kmer_iterator{_view, next};
        }

        reference operator*()
        {
            return _kmer;
        }

        pointer operator->()
        {
            return &_kmer;
        }

    private:
        const branch_entry_view& _view;
        phylo_kmer _kmer;
        size_t _start_pos;
        size_t _last_letter_idx;
    };

public:
    using iterator = phylo_kmer_iterator<false>;
    using const_iterator = const phylo_kmer_iterator<true>;
    using reference = iterator::reference;
    using const_reference = const_iterator::reference;

    branch_entry_view(const branch_entry& entry, size_t begin, size_t end);
    branch_entry_view(const branch_entry_view&) = default;
    branch_entry_view(branch_entry_view&&) = delete;
    branch_entry_view& operator=(const branch_entry_view&) = delete;
    branch_entry_view&& operator=(branch_entry_view&&) = delete;
    ~branch_entry_view() = default;

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    size_t kmer_size() const;
    size_t alphabet_size() const;

    const proba_pair& at(size_t position, size_t variant) const;

private:
    const branch_entry& _entry;
    const size_t _begin;
    const size_t _end;
};

#endif //RAPPAS_CPP_BRANCH_ENTRY_H
