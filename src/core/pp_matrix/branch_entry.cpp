#include "branch_entry.h"

branch_entry_view::branch_entry_view(const branch_entry* entry, size_t start, size_t end) noexcept
    : _entry{ entry }
    , _start{ start }
    , _end{ end }
{}

branch_entry_view::branch_entry_view(const branch_entry_view& other) noexcept
    : branch_entry_view{ other._entry, other._start, other._end}
{}


branch_entry_view::iterator branch_entry_view::begin()
{
    return phylo_kmer_iterator<false>{ this, _end - _start, _start, first_phylo_kmer() };
}

branch_entry_view::iterator branch_entry_view::end()
{
    return phylo_kmer_iterator<false>{ this, 0, 0, {}};
}

branch_entry_view::const_iterator branch_entry_view::begin() const
{
    return phylo_kmer_iterator<true>{ this, _end - _start, _start, first_phylo_kmer() };
}

branch_entry_view::const_iterator branch_entry_view::end() const
{
    return phylo_kmer_iterator<true>{ this, 0, 0, {}};
}

const proba_pair& branch_entry_view::at(size_t position, size_t variant) const
{
    return _entry->at(position, variant);
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

size_t branch_entry_view::get_start_pos() const
{
    return _start;
}

size_t branch_entry_view::get_end_pos() const
{
    return _end;
}

phylo_kmer_iterator<false>::stack_type branch_entry_view::first_phylo_kmer() const
{
    const auto kmer_size = _end - _start;
    const auto alphabet_size = _entry->get_alphabet_size();

    auto stack = phylo_kmer_iterator<true>::stack_type{};
    phylo_kmer kmer { 0, 0.0f };
    for (size_t i = 0; i < kmer_size; ++i)
    {
        const size_t next_move = 0;
        const auto& proba_pair = _entry->at(i, next_move);
        kmer.value = left_shift(kmer.value, kmer_size, alphabet_size) | proba_pair.index;
        kmer.score += proba_pair.score;
        stack.emplace_back(next_move, kmer);
    }
    return stack;
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

branch_entry::branch_entry(const seq_traits& traits)
    : _traits{ traits }
{}

branch_entry::branch_entry(branch_id _id, std::vector<row>&& rows, const seq_traits& traits)
    : _branch_id{ _id }
    , _rows{ std::move(rows) }
    , _traits{ traits }
{}

branch_entry::iterator branch_entry::begin(size_t kmer_size)
{
    return { {this, 0, kmer_size } };
}

branch_entry::iterator branch_entry::end()
{
    return { {this, 0, 0} };
}

branch_entry::const_iterator branch_entry::begin(size_t kmer_size) const
{
    return { {this, 0, kmer_size } };
}

branch_entry::const_iterator branch_entry::end() const
{
    return { {this, 0, 0} };
}

void branch_entry::push_back(row&& r)
{
    _rows.push_back(std::move(r));
}

size_t branch_entry::get_alignment_size() const
{
    return _rows.size();
}

size_t branch_entry::get_alphabet_size() const
{
    return std::begin(_rows)->size();
}

branch_id branch_entry::get_branch_id() const
{
    return _branch_id;
}

const seq_traits& branch_entry::traits() const
{
    return _traits;
}

const proba_pair& branch_entry::at(size_t position, size_t variant) const
{
    return _rows[position][variant];
}

bool operator==(const branch_entry& lhs, const branch_entry& rhs)
{
    return lhs.get_branch_id() == rhs.get_branch_id();
}