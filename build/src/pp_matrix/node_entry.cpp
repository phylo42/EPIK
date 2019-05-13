#include "node_entry.h"

using namespace rappas;

node_entry::node_entry(branch_type _id, vector_type&& rows)
    : _branch_label{ _id }
    , _rows{ std::move(rows) }
{}

node_entry::const_iterator node_entry::begin(uint32_t kmer_size) const
{
    return { { this, 0, kmer_size } };
}

node_entry::const_iterator node_entry::end() const
{
    return { { this, 0, 0 } };
}

void node_entry::push_back(row_type&& row)
{
    _rows.push_back(row);
}

size_t node_entry::get_alignment_size() const
{
    return _rows.size();
}

branch_type node_entry::get_label() const
{
    return _branch_label;
}

const proba_pair& node_entry::at(size_t position, size_t variant) const
{
    return _rows[position][variant];
}

bool operator==(const node_entry& lhs, const node_entry& rhs)
{
    return lhs.get_label() == rhs.get_label();
}

view_iterator::view_iterator(node_entry_view view) noexcept
    : _view{ view }
{}

view_iterator& view_iterator::operator++()
{
    auto entry = _view.get_entry();
    if (_view.get_end_pos() < entry->get_alignment_size())
    {
        _view = { node_entry_view{ entry, _view.get_start_pos() + 1, _view.get_end_pos() + 1 } };
    }
    else
    {
        /// WARNING: the same code in node_entry::end
        _view = { node_entry_view{ entry, 0, 0 } };
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