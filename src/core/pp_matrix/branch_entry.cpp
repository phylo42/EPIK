#include "branch_entry.h"

branch_entry::branch_entry(const seq_traits& traits)
    : branch_entry{{}, traits}
{}

branch_entry::branch_entry(std::vector<row>&& rows, const seq_traits& traits)
    : _rows{std::move(rows)}
    , _traits{traits}
{}

void branch_entry::push_back(row&& r)
{
    _rows.push_back(std::move(r));
}

size_t branch_entry::alignment_size() const
{
    return _rows.size();
}

size_t branch_entry::alphabet_size() const
{
    return std::begin(_rows)->size();
}

const proba_pair& branch_entry::at(size_t position, size_t variant) const
{
    return _rows[position][variant];
}
branch_entry_view::branch_entry_view(const branch_entry& entry, size_t begin, size_t end)
    : _entry(entry)
    , _begin(begin)
    , _end(end)
{}

const proba_pair& branch_entry_view::at(size_t position, size_t variant) const
{
    return _entry.at(position, variant);
}

size_t branch_entry_view::kmer_size() const
{
    return _end - _begin;
}

size_t branch_entry_view::alphabet_size() const
{
    return _entry.alphabet_size();
}