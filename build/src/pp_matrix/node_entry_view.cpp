#include "node_entry.h"
#include "node_entry_view.h"
#include <cmath>
#include <algorithm>
#include <xpas/seq.h>

using namespace xpas;
using namespace rappas::impl;

dac_kmer_iterator make_dac_end_iterator()
{
    return { nullptr, 0, 0, 0 };
}

bool kmer_score_comparator(const unpositioned_phylo_kmer& k1, const unpositioned_phylo_kmer& k2)
{
    return k1.score > k2.score;
}

dac_kmer_iterator::dac_kmer_iterator(const node_entry* entry, size_t kmer_size, xpas::phylo_kmer::score_type threshold,
                                     xpas::phylo_kmer::pos_type start_pos) noexcept
    : _entry{ entry }, _kmer_size{ kmer_size }, _left_part_size{ 0 }, _start_pos{ start_pos }, _threshold{ threshold }
{
    const auto halfsize = size_t{ kmer_size / 2 };
    _left_part_size = (halfsize >= 1)
        ? size_t{ kmer_size / 2 }
        : kmer_size;

    if (kmer_size == 1)
    {
        for (size_t i = 0; i < seq_traits::alphabet_size; ++i)
        {
            const auto& letter = _entry->at(_start_pos, i);
            _left_halfmers.push_back(make_phylo_kmer<unpositioned_phylo_kmer>(letter.index, letter.score, 0));
        }
        _left_halfmer_it = _left_halfmers.begin();
        _current = _next_phylokmer();
    }
    /// a workaround for end()
    else if (_left_part_size > 0)
    {
        // left part
        {
            auto it = (_left_part_size == 0)
                      ? make_dac_end_iterator()
                      : dac_kmer_iterator(entry, _left_part_size, threshold, start_pos);
            const auto end = make_dac_end_iterator();
            for (; it != end; ++it)
            {
                _left_halfmers.push_back(*it);
            }
            _left_halfmer_it = _left_halfmers.begin();
        }

        // right part
        {
            auto it = (_left_part_size < kmer_size)
                      ? dac_kmer_iterator(entry, kmer_size - _left_part_size, threshold, start_pos + _left_part_size)
                      : make_dac_end_iterator();
            const auto end = make_dac_end_iterator();
            for (; it != end; ++it)
            {
                _right_halfmers.push_back(*it);
            }
            std::sort(_right_halfmers.begin(), _right_halfmers.end(), kmer_score_comparator);
            _right_halfmer_it = _right_halfmers.begin();
            _select_right_halfmers_bound();
        }

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
        _left_halfmers = std::move(rhs._left_halfmers);
        _left_halfmer_it = rhs._left_halfmer_it;
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

unpositioned_phylo_kmer dac_kmer_iterator::_next_phylokmer()
{
    if (_left_halfmer_it != _left_halfmers.end())
    {
        if (_kmer_size > 1)
        {
            while ((_right_halfmer_it == _last_right_halfmer_it) && (_left_halfmer_it != _left_halfmers.end()))
            {
                ++_left_halfmer_it;
                _right_halfmer_it = _right_halfmers.begin();
                _select_right_halfmers_bound();
            }

            if (_left_halfmer_it != _left_halfmers.end())
            {
                const auto left_halfmer = *_left_halfmer_it;
                const auto right_halfmer = *_right_halfmer_it;
                const auto full_key =
                    (left_halfmer.key << ((_kmer_size - _left_part_size) * xpas::bit_length<xpas::seq_type>()))
                    | right_halfmer.key;
                const auto full_score = left_halfmer.score + right_halfmer.score;
                ++_right_halfmer_it;
                return make_phylo_kmer<unpositioned_phylo_kmer>(full_key, full_score, 0);
            }
        }
        else
        {
            const auto kmer = *_left_halfmer_it;
            ++_left_halfmer_it;
            return kmer;
        }
    }

    *this = make_dac_end_iterator();
    return {};
}

void dac_kmer_iterator::_select_right_halfmers_bound()
{
    const auto residual_threshold =  _threshold - _left_halfmer_it->score;
    _last_right_halfmer_it = ::std::lower_bound(_right_halfmers.begin(), _right_halfmers.end(),
        make_phylo_kmer<unpositioned_phylo_kmer>(0, residual_threshold, 0), kmer_score_comparator);
}

node_entry_view::node_entry_view(const node_entry* entry, phylo_kmer::score_type threshold,
                                 phylo_kmer::pos_type start, phylo_kmer::pos_type end) noexcept
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
    return make_dac_end_iterator();
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

xpas::phylo_kmer::pos_type node_entry_view::get_start_pos() const noexcept
{
    return _start;
}

xpas::phylo_kmer::pos_type node_entry_view::get_end_pos() const noexcept
{
    return _end;
}

xpas::phylo_kmer::score_type node_entry_view::get_threshold() const noexcept
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
