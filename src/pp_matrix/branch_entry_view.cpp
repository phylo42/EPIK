#include "branch_entry.h"
#include "branch_entry_view.h"
#include <cmath>

phylo_kmer_iterator::phylo_kmer_iterator(const branch_entry* entry, size_t kmer_size,
                                         core::phylo_kmer::pos_type start_pos, stack_type stack) noexcept
    : _entry{ entry }
    , _kmer_size{ kmer_size }
    , _start_pos{ start_pos }
    , _threshold{ core::score_threshold(kmer_size) }
    , _stack{ std::move(stack) }
{}

bool phylo_kmer_iterator::operator==(const phylo_kmer_iterator& rhs) const noexcept
{
    if (_stack.empty())
    {
        return rhs._stack.empty();
    }
    if (rhs._stack.empty())
    {
        return false;
    }
    return _stack.back().mmer == rhs._stack.back().mmer;
}

bool phylo_kmer_iterator::operator!=(const phylo_kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

phylo_kmer_iterator& phylo_kmer_iterator::operator++()
{
    _current = next_phylokmer();
    return *this;
}

phylo_kmer_iterator::reference phylo_kmer_iterator::operator*()
{
    return _current.mmer;
}

phylo_kmer_iterator::pointer phylo_kmer_iterator::operator->()
{
    return &(_current.mmer);
}

void phylo_kmer_iterator::next_index()
{
    const auto last_mmer = _stack.back();
    _stack.pop_back();
    if (last_mmer.last_index < core::seq_traits::alphabet_size - 1)
    {
        const auto new_letter_position = last_mmer.last_position;
        const auto new_letter_index = last_mmer.last_index + 1;
        const auto last_letter_score = _entry->at(_start_pos + last_mmer.last_position, last_mmer.last_index).score;
        const auto& new_letter = _entry->at(_start_pos + new_letter_position, new_letter_index);
        const auto new_mmer_key = (last_mmer.mmer.key & core::rightest_symbol_mask<core::seq_type>()) | new_letter.index;
        const auto new_mmer_score = last_mmer.mmer.score - last_letter_score + new_letter.score;
        const auto new_mmer = phylo_mmer{ {new_mmer_key, new_mmer_score}, new_letter_position, new_letter_index, false };
        _stack.push_back(new_mmer);
    }
    else
    {
        if (!_stack.empty())
        {
            _stack.back().visited = true;
        }
    }
}

void phylo_kmer_iterator::next_position()
{
    const auto& last_mmer = _stack.back();
    const auto new_letter_position = last_mmer.last_position + 1;
    const auto new_letter_index = 0;
    const auto& new_letter = _entry->at(_start_pos + new_letter_position, new_letter_index);
    const auto new_mmer_key = (last_mmer.mmer.key << core::bit_length<core::seq_type>()) | new_letter.index;
    const auto new_mmer_score = last_mmer.mmer.score + new_letter.score;
    const auto new_mmer = phylo_mmer{ {new_mmer_key, new_mmer_score}, new_letter_position, new_letter_index, false };
    _stack.push_back(new_mmer);
}

phylo_kmer_iterator::phylo_mmer phylo_kmer_iterator::next_phylokmer()
{
    while (!_stack.empty())
    {
        /// cut a subtree (and go bottom-up)
        if (_stack.back().mmer.score <= _threshold)
        {
            _stack.pop_back();
            if (!_stack.empty())
            {
                _stack.back().visited = true;
            }
        }
        /// go top-down the tree
        else if (_stack.size() < _kmer_size)
        {
            if (_stack.back().visited)
            {
                next_index();
            }
            else
            {
                next_position();
            }
        }
        /// we are at the leaf level, test a k-mer at next iteration or return
        else
        {
            const auto last_mmer = _stack.back();
            if (last_mmer.visited)
            {
                next_index();
                return last_mmer;
            }
            else
            {
                _stack.back().visited = true;
            }
        }
    }
    return {};
}

branch_entry_view::branch_entry_view(const branch_entry* entry, core::phylo_kmer::pos_type start, core::phylo_kmer::pos_type end) noexcept
    : _entry{ entry }
    , _start{ start }
    , _end{ end }
{}

branch_entry_view::branch_entry_view(const branch_entry_view& other) noexcept
    : branch_entry_view{ other._entry, other._start, other._end}
{}

branch_entry_view::const_iterator branch_entry_view::begin() const
{
    const auto& first_cell = _entry->at(_start, 0);
    phylo_kmer_iterator::stack_type stack;
    stack.push_back(phylo_kmer_iterator::phylo_mmer{ { first_cell.index, first_cell.score }, 0, 0, false });
    auto it = phylo_kmer_iterator{ _entry, _end - _start, _start, stack };
    ++it;
    return it;
}

branch_entry_view::const_iterator branch_entry_view::end() const
{
    return phylo_kmer_iterator{ _entry, 0, 0, {}};
}

branch_entry_view& branch_entry_view::operator=(branch_entry_view&& other) noexcept
{
    if (*this != other)
    {
        _entry = other._entry;
        _start = other._start;
        _end = other._end;
        other._entry = nullptr;
        other._start = 0;
        other._end = 0;
    }
    return *this;
}

const branch_entry* branch_entry_view::get_entry() const
{
    return _entry;
}

core::phylo_kmer::pos_type branch_entry_view::get_start_pos() const
{
    return _start;
}

core::phylo_kmer::pos_type branch_entry_view::get_end_pos() const
{
    return _end;
}

bool operator==(const branch_entry_view& a, const branch_entry_view& b) noexcept
{
    return (a.get_start_pos() == b.get_start_pos()) &&
           (a.get_end_pos() == b.get_end_pos()) &&
           (a.get_entry() == b.get_entry());
}

bool operator!=(const branch_entry_view& a, const branch_entry_view& b) noexcept
{
    return !(a == b);
}
