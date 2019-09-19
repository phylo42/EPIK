#include "node_entry.h"
#include "node_entry_view.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <core/seq.h>

using namespace core;
using namespace rappas::impl;


bnb_kmer_iterator make_bnb_end_iterator()
{
    return bnb_kmer_iterator();
}

bnb_kmer_iterator make_bnb_begin_iterator(const node_entry* entry, phylo_kmer::pos_type start, size_t kmer_size,
                                          phylo_kmer::score_type threshold)
{
    bnb_kmer_iterator::stack_type stack;
    stack.reserve(kmer_size + 1);

    /// Push a fake mmer to the bottom of the stack. We use it to reduce the number
    /// of if statements in the operator++, and therefore to reduce the number of branch
    /// mispredictions. It is also necessary to start branch-and-bound with every value
    /// of the first column without adding new if statements.
    stack.push_back(phylo_mmer{ { 0, 0.0 }, -1, 0, 0 });

    /// Calculate the first k-mer
    phylo_kmer::key_type kmer_key = 0;
    phylo_kmer::score_type kmer_score = 0.0;
    for (size_t i = 0; i < kmer_size; ++i)
    {
        const auto& ith_letter = entry->at(start + i, 0);
        kmer_key = (kmer_key << bit_length<seq_type>()) | ith_letter.index;
        kmer_score += ith_letter.score;

        stack.push_back(phylo_mmer{ { kmer_key, kmer_score }, phylo_kmer::pos_type(i), 0, 0 });
    }

    /// There is no point to iterate over the window if the first k-mer is already not good enough
    if (kmer_score < threshold)
    {
        return make_bnb_end_iterator();
    }
    else
    {
        return bnb_kmer_iterator{ entry, kmer_size, threshold, start, std::move(stack) };
    }
}

bnb_kmer_iterator::bnb_kmer_iterator() noexcept
    : _entry{ nullptr }, _kmer_size{ 0 }, _start_pos{ 0 }, _threshold{ 0.0 }, _stack{ }
{}

bnb_kmer_iterator::bnb_kmer_iterator(const node_entry* entry, size_t kmer_size, phylo_kmer::score_type threshold,
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

bnb_kmer_iterator& bnb_kmer_iterator::operator=(bnb_kmer_iterator&& rhs) noexcept
{
    if (*this != rhs)
    {
        _entry = rhs._entry;
        _kmer_size = rhs._kmer_size;
        _start_pos = rhs._start_pos;
        _threshold = rhs._threshold;
        _stack = std::move(rhs._stack);
        _current = rhs._current;
    }
    return *this;
}

bool bnb_kmer_iterator::operator==(const bnb_kmer_iterator& rhs) const noexcept
{
    if (_entry != rhs._entry)
    {
        return false;
    }
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

bool bnb_kmer_iterator::operator!=(const bnb_kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

bnb_kmer_iterator& bnb_kmer_iterator::operator++()
{
    _current = next_phylokmer();
    return *this;
}

bnb_kmer_iterator::reference bnb_kmer_iterator::operator*() const noexcept
{
    return _current.mmer;
}

bnb_kmer_iterator::pointer bnb_kmer_iterator::operator->() const noexcept
{
    return &(_current.mmer);
}

phylo_mmer bnb_kmer_iterator::next_phylokmer()
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

    *this = make_bnb_end_iterator();
    return {};
}

dac_kmer_iterator make_daq_end_iterator()
{
    return { nullptr, 0, 0, 0 };
}

bool kmer_score_comparator(const phylo_kmer& k1, const phylo_kmer& k2)
{
    return k1.score > k2.score;
}

dac_kmer_iterator::dac_kmer_iterator(const node_entry* entry, size_t kmer_size, core::phylo_kmer::score_type threshold,
                                     core::phylo_kmer::pos_type start_pos) noexcept
    : _entry{ entry }, _kmer_size{ kmer_size }, _left_part_size{ 0 }, _start_pos{ start_pos }, _threshold{ threshold }
{
    const auto halfsize = size_t{ kmer_size / 2 };
    _left_part_size = (halfsize >= 1)
        ? size_t{ kmer_size / 2 }
        : kmer_size;

    _left_iterator = make_bnb_begin_iterator(entry, start_pos, _left_part_size, threshold);

    /// a workaround for end()
    if (_left_part_size > 0)
    {
        _right_halfmers.reserve(size_t(std::pow(seq_traits::alphabet_size, _kmer_size - _left_part_size)) / 2);
        auto it = (_left_part_size < kmer_size)
            ? make_bnb_begin_iterator(entry, start_pos + _left_part_size, kmer_size - _left_part_size, threshold)
            : make_bnb_end_iterator();
        const auto end = make_bnb_end_iterator();
        for (; it != end; ++it)
        {
            _right_halfmers.push_back(*it);
        }
        std::sort(_right_halfmers.begin(), _right_halfmers.end(), kmer_score_comparator);
        _right_halfmer_it = _right_halfmers.begin();
        _select_right_halfmers_bound();

        _current = _next_phylokmer();
    }
}

dac_kmer_iterator& dac_kmer_iterator::operator=(dac_kmer_iterator&& rhs) noexcept
{
    if (*this != rhs)
    {
        _entry = rhs._entry;
        _kmer_size = rhs._kmer_size;
        _left_part_size = rhs._left_part_size;
        _start_pos = rhs._start_pos;
        _threshold = rhs._threshold;
        _current = rhs._current;
        _left_iterator = std::move(rhs._left_iterator);
        _right_halfmers = std::move(rhs._right_halfmers);
        _right_halfmer_it = rhs._right_halfmer_it;
        _last_right_halfmer_it = rhs._last_right_halfmer_it;
    }
    return *this;
}

bool dac_kmer_iterator::operator==(const dac_kmer_iterator& rhs) const noexcept
{
    return _entry == rhs._entry && _start_pos == rhs._start_pos && _kmer_size == rhs._kmer_size;
}

bool dac_kmer_iterator::operator!=(const dac_kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

dac_kmer_iterator& dac_kmer_iterator::operator++()
{
    _current = _next_phylokmer();
    return *this;
}

dac_kmer_iterator::reference dac_kmer_iterator::operator*() const noexcept
{
    return _current;
}

dac_kmer_iterator::pointer dac_kmer_iterator::operator->()  const noexcept
{
    return &_current;
}

phylo_kmer dac_kmer_iterator::_next_phylokmer()
{
    while (_right_halfmer_it == _last_right_halfmer_it)
    {
        ++_left_iterator;
        _right_halfmer_it = _right_halfmers.begin();
        _select_right_halfmers_bound();
    }

    if (_left_iterator != make_bnb_end_iterator())
    {
        const auto left_halfmer = *_left_iterator;
        const auto right_halfmer = *_right_halfmer_it;
        const auto full_key = (left_halfmer.key << ((_kmer_size - _left_part_size) * core::bit_length<core::seq_type>()))
                              | right_halfmer.key;
        const auto full_score = left_halfmer.score + right_halfmer.score;
        ++_right_halfmer_it;
        return { full_key, full_score };
    }
    else
    {
        *this = make_daq_end_iterator();
        return {};
    }
}

void dac_kmer_iterator::_select_right_halfmers_bound()
{
    const auto residual_threshold =  _threshold - _left_iterator->score;
    _last_right_halfmer_it = ::std::lower_bound(_right_halfmers.begin(), _right_halfmers.end(),
        phylo_kmer{ 0, residual_threshold }, kmer_score_comparator);
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
    return { _entry, kmer_size, _threshold, _start };
}

node_entry_view::const_iterator node_entry_view::end() const noexcept
{
    return make_daq_end_iterator();
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
