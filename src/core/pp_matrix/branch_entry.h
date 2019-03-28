#ifndef RAPPAS_CPP_BRANCH_ENTRY_H
#define RAPPAS_CPP_BRANCH_ENTRY_H

#include <vector>
#include <core/seq_traits.h>
#include <utils/meta.h>
#include "row.h"

class branch_entry;
class branch_entry_view;

/// \brief A forward access (non-)const iterator for phylo_kmer pairs [kmer value, score]. Iterates over
/// a fixed branch_entry_view of size K.
template <bool IsConst>
class phylo_kmer_iterator
{
public:
    using iterator_category = std::forward_iterator_tag;
    using reference = typename choose<IsConst, const phylo_kmer&, phylo_kmer&>::type;
    using pointer = typename choose<IsConst, const phylo_kmer*, phylo_kmer*>::type;
    using stack_type = std::vector<std::pair<size_t, phylo_kmer>>;

    phylo_kmer_iterator(const branch_entry_view* view, size_t kmer_size,
                        size_t start_pos, stack_type stack) noexcept
        : _view{ view }
        , _kmer_size{ kmer_size }
        , _start_pos{ start_pos }
        , _stack{ std::move(stack) }
    {}
    phylo_kmer_iterator(const phylo_kmer_iterator&) = delete;
    phylo_kmer_iterator(phylo_kmer_iterator&&) noexcept = default;
    phylo_kmer_iterator& operator=(const phylo_kmer_iterator& rhs) = default;
    phylo_kmer_iterator& operator=(phylo_kmer_iterator&&) noexcept = default;
    ~phylo_kmer_iterator() noexcept = default;

    bool operator==(const phylo_kmer_iterator& rhs) const noexcept
    {
        if (_stack.empty())
        {
            return rhs._stack.empty();
        }
        if (rhs._stack.empty())
        {
            return false;
        }
        return _stack.back() == rhs._stack.back();
    }

    bool operator!=(const phylo_kmer_iterator& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    phylo_kmer_iterator& operator++();

    reference operator*()
    {
        return _stack.back().second;
    }

    pointer operator->()
    {
        return &(_stack.back().second);
    }

private:
    const branch_entry_view* _view;
    size_t _kmer_size;
    size_t _start_pos;
    stack_type _stack;
};

/// \brief A lightweight view of branch_entry. Implements a "window" of size K over a branch_entry.
class branch_entry_view final
{
public:
    using iterator = phylo_kmer_iterator<false>;
    using const_iterator = const phylo_kmer_iterator<true>;
    using reference = iterator::reference;
    using const_reference = const_iterator::reference;

    branch_entry_view(const branch_entry* entry, size_t start, size_t end) noexcept;
    branch_entry_view(const branch_entry_view& other) noexcept;
    branch_entry_view(branch_entry_view&&) = delete;
    branch_entry_view& operator=(const branch_entry_view&) = delete;
    branch_entry_view& operator=(branch_entry_view&& other) noexcept;
    ~branch_entry_view() noexcept = default;

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

    const proba_pair& at(size_t position, size_t variant) const;

    const branch_entry* get_entry() const;
    size_t get_start_pos() const;
    size_t get_end_pos() const;

private:
    phylo_kmer_iterator<false>::stack_type first_phylo_kmer() const;

    const branch_entry* _entry;
    size_t _start;
    size_t _end;
};

bool operator==(const branch_entry_view& a, const branch_entry_view& b) noexcept;
bool operator!=(const branch_entry_view& a, const branch_entry_view& b) noexcept;


template <bool IsConst>
class view_iterator;

/// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
class branch_entry final
{
public:
    using iterator = view_iterator<false>;
    using const_iterator = view_iterator<true>;

    explicit branch_entry(const seq_traits& traits = dna_seq_traits);
    branch_entry(branch_id _id, std::vector<row>&& rows, const seq_traits& traits);
    branch_entry(const branch_entry&) = delete;
    branch_entry(branch_entry&&) = default;
    branch_entry& operator=(const branch_entry&) = delete;
    branch_entry& operator=(branch_entry&&) = default;
    ~branch_entry() noexcept = default;

    iterator begin(size_t kmer_size);
    iterator end();
    const_iterator begin(size_t kmer_size) const;
    const_iterator end() const;

    void push_back(row&& r);

    size_t get_alignment_size() const;
    size_t get_alphabet_size() const;
    branch_id get_branch_id() const;

    const seq_traits& traits() const;

    const proba_pair& at(size_t position, size_t variant) const;

private:
    branch_id _branch_id;
    std::vector<row> _rows;
    seq_traits _traits;
};

bool operator==(const branch_entry& lhs, const branch_entry& rhs);

template <bool IsConst>
class view_iterator
{
public:
    typedef std::forward_iterator_tag iterator_category;
    typedef typename choose<IsConst, const branch_entry_view&, branch_entry_view&>::type reference;
    typedef typename choose<IsConst, const branch_entry_view*, branch_entry_view*>::type pointer;

    view_iterator(branch_entry_view view) : _view{ view } {}
    view_iterator(const view_iterator& view) = delete;
    view_iterator(view_iterator&& view) = delete;
    view_iterator& operator=(const view_iterator&) = delete;
    view_iterator& operator=(view_iterator&&) = delete;
    ~view_iterator() = default;

    view_iterator& operator++()
    {
        auto entry = _view.get_entry();
        if (_view.get_end_pos() < entry->get_alignment_size() - 1)
        {
            _view = { branch_entry_view{ entry, _view.get_start_pos() + 1, _view.get_end_pos() + 1 } };
        }
        else
        {
            /// WARNING: the same code in branch_entry::end
            _view = { branch_entry_view{ entry, 0, 0 } };
        }
        return *this;
    }

    bool operator==(const view_iterator& rhs) const noexcept
    {
        return _view == rhs._view;
    }

    bool operator!=(const view_iterator& rhs) const noexcept
    {
        return !(*this == rhs);
    }

    reference operator*()
    {
        return _view;
    }

    pointer operator->()
    {
        return &_view;
    }

private:
    branch_entry_view _view;
};

template<bool IsConst>
phylo_kmer_iterator<IsConst>& phylo_kmer_iterator<IsConst>::operator++()
{
    const auto entry = _view->get_entry();
    const auto alphabet_size = entry->get_alphabet_size();

    /// TODO: measure if copying phylo_kmer struct is effective here
    auto& last_pair = _stack.back();
    auto last_move = last_pair.first;
    auto last_letter_row = _start_pos + _stack.size() - 1;
    auto last_letter = _view->at(last_letter_row, last_move);
    auto& last_kmer = last_pair.second;
    auto next_kmer_value = last_kmer.value;
    auto next_kmer_score = last_kmer.score;
    auto next_move = last_move + 1;

    /// go bottom-up and pop out m-mers (m <= k)
    do
    {
        /// if we discovered the last k-mer
        if (_stack.empty())
        {
            *this = _view->end();
            return *this;
        }

        _stack.pop_back();

        /// if we have k-1-mer, remember which letter we need to add, otherwise we add first letter
        if (_stack.size() < _kmer_size)
        {
            next_move = last_move + 1;
        }
        else
        {
            next_move = 0;
        }

        if (_stack.empty())
        {
            last_letter_row = _start_pos;
        }
        else
        {
            last_pair = _stack.back();
            last_move = last_pair.first;
            last_kmer = last_pair.second;

            last_letter_row = _start_pos + _stack.size() - 1;
            last_letter = _view->at(last_letter_row, last_move);

            next_kmer_value = last_kmer.value;
            next_kmer_score = last_kmer.score;
        }
    } while (next_move >= alphabet_size);

    /// go top-down and push undiscovered m-mers (m <= k)
    const auto kmask = mask(_kmer_size, alphabet_size);
    while (_stack.size() < _kmer_size)
    {
        /// calculate new k-mer value and score
        auto new_letter = _view->at(last_letter_row + 1, next_move);

        next_kmer_value = (left_shift(next_kmer_value, _kmer_size, alphabet_size) & kmask) | new_letter.index;
        next_kmer_score += new_letter.score;
        _stack.emplace_back(next_move, phylo_kmer{ next_kmer_value, next_kmer_score });

        ++last_letter_row;
        next_move = 0;
    }

    return *this;
}
#endif //RAPPAS_CPP_BRANCH_ENTRY_H
