#include "branch_entry.h"

branch_entry::branch_entry(branch_id _id, vector_type&& rows)
    : _branch_label{ _id }
    , _rows{ std::move(rows) }
{}

branch_entry::const_iterator branch_entry::begin(uint32_t kmer_size) const
{
    return { { this, 0, kmer_size } };
}

branch_entry::const_iterator branch_entry::end() const
{
    return { { this, 0, 0 } };
}

void branch_entry::push_back(row&& r)
{
    _rows.push_back(r);
}

size_t branch_entry::get_alignment_size() const
{
    return _rows.size();
}

branch_id branch_entry::get_branch_label() const
{
    return _branch_label;
}

const proba_pair& branch_entry::at(size_t position, size_t variant) const
{
    return _rows[position][variant];
}

bool operator==(const branch_entry& lhs, const branch_entry& rhs)
{
    return lhs.get_branch_label() == rhs.get_branch_label();
}

view_iterator::view_iterator(branch_entry_view view) noexcept
    : _view{ view }
{}

view_iterator& view_iterator::operator++()
{
    auto entry = _view.get_entry();
    if (_view.get_end_pos() < entry->get_alignment_size())
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

bool view_iterator::operator==(const view_iterator& rhs) const noexcept
{
    return _view == rhs._view;
}

bool view_iterator::operator!=(const view_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

view_iterator::const_reference view_iterator::operator*() const noexcept
{
    return _view;
}

view_iterator::const_pointer view_iterator::operator->() const noexcept
{
    return &_view;
}