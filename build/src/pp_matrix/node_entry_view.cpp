#include "node_entry.h"
#include "node_entry_view.h"
#include <cmath>
#include <iostream>

using namespace core;

phylo_kmer_iterator::phylo_kmer_iterator(const node_entry* entry, size_t kmer_size, phylo_kmer::score_type threshold,
                                         phylo_kmer::pos_type start_pos, stack_type&& stack) noexcept
    : _entry{ entry }
    , _kmer_size{ kmer_size }
    , _start_pos{ start_pos }
    , _threshold{ threshold }
    , _stack{ std::move(stack) }
{
    /// stack can be empty for the end() method
    if (!_stack.empty())
    {
        _current = _stack.back();
    }
}

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

phylo_kmer_iterator::phylo_mmer phylo_kmer_iterator::next_phylokmer()
{
    {
        const auto next_index = _stack.back().last_index + 1;
        _stack.pop_back();
        _stack.back().next_index = next_index;
    }

    /// until there is only the fake k-mer
    while (!(_stack.size() == 1 && _stack.back().next_index == core::seq_traits::alphabet_size))
    {
        const auto top_mmer = _stack.back();

        /// if we found a m-mer with low score, we cut the whole subtree
        if (top_mmer.mmer.score < _threshold)
        {
            _stack.pop_back();
            _stack.back().next_index = core::seq_traits::alphabet_size;
        }
        else
        {
            /// |top_mmer| < k, go to the next column
            /// + 1 because of the fake mmer
            if (_stack.size() < _kmer_size + 1)
            {
                /// go to the next letter in the next column
                if (top_mmer.next_index < core::seq_traits::alphabet_size)
                {
                    const auto new_letter_position = top_mmer.last_position + 1;
                    const auto& new_letter = _entry->at(_start_pos + new_letter_position, top_mmer.next_index);
                    const auto new_mmer_key = (top_mmer.mmer.key << core::bit_length<core::seq_type>()) | new_letter.index;
                    const auto new_mmer_score = top_mmer.mmer.score + new_letter.score;
                    _stack.push_back(phylo_mmer{ { new_mmer_key, new_mmer_score }, new_letter_position, top_mmer.next_index, 0 });
                }
                /// the next column is over, get back
                else
                {
                    _stack.pop_back();
                    _stack.back().next_index = top_mmer.last_index + 1;
                }
            }
            /// we have a good k-mer
            else
            {
                return top_mmer;
            }
        }
    }
    /// to clear the stack literally means
    /// "*this = node_entry_view::end()"
    _stack.clear();
    return {};
}

node_entry_view::node_entry_view(const node_entry* entry, core::phylo_kmer::score_type threshold,
    core::phylo_kmer::pos_type start, core::phylo_kmer::pos_type end) noexcept
    : _entry{ entry }, _threshold{ threshold }, _start{ start }, _end{ end }
{}

node_entry_view::node_entry_view(const node_entry_view& other) noexcept
    : node_entry_view{ other._entry, other._threshold, other._start, other._end}
{}

node_entry_view::const_iterator node_entry_view::begin() const
{
    const auto kmer_size = size_t{ (size_t)_end - _start };
    phylo_kmer_iterator::stack_type stack;

    /// Push a fake mmer to the bottom of the stack. We use it to reduce the number
    /// of if statements in the operator++, and therefore to reduce the number of branch
    /// mispredictions. It is also necessary to start branch-and-bound with every value
    /// of the first column without adding new if statements.
    stack.push_back(phylo_kmer_iterator::phylo_mmer{ { 0, 0.0 }, -1, 0, 0 });

    /// Calculate the first k-mer
    phylo_kmer::key_type kmer_key = 0;
    phylo_kmer::score_type kmer_score = 0.0;
    for (size_t i = 0; i < kmer_size; ++i)
    {
        const auto& ith_letter = _entry->at(_start + i, 0);
        kmer_key = (kmer_key << core::bit_length<core::seq_type>()) | ith_letter.index;
        kmer_score += ith_letter.score;

        stack.push_back(phylo_kmer_iterator::phylo_mmer{ { kmer_key, kmer_score }, core::phylo_kmer::pos_type(i), 0, 0 });
    }

    /// There is no point to iterate over the window if the first k-mer is already not good enough
    if (kmer_score < _threshold)
    {
        return end();
    }
    else
    {
        return phylo_kmer_iterator{ _entry, kmer_size, _threshold, _start, std::move(stack) };
    }
}

node_entry_view::const_iterator node_entry_view::end() const noexcept
{
    return phylo_kmer_iterator{ _entry, 0, _threshold, 0, {}};
}

node_entry_view& node_entry_view::operator=(node_entry_view&& other) noexcept
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

const node_entry* node_entry_view::get_entry() const noexcept
{
    return _entry;
}

core::phylo_kmer::pos_type node_entry_view::get_start_pos() const noexcept
{
    return _start;
}

core::phylo_kmer::pos_type node_entry_view::get_end_pos() const noexcept
{
    return _end;
}

core::phylo_kmer::score_type node_entry_view::get_threshold() const noexcept
{
    return _threshold;
}

bool operator==(const node_entry_view& a, const node_entry_view& b) noexcept
{
    return (a.get_start_pos() == b.get_start_pos()) &&
           (a.get_end_pos() == b.get_end_pos()) &&
           (a.get_entry() == b.get_entry());
}

bool operator!=(const node_entry_view& a, const node_entry_view& b) noexcept
{
    return !(a == b);
}
