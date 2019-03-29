#include "branch_entry.h"
#include "branch_entry_view.h"

phylo_kmer_iterator::phylo_kmer_iterator(const branch_entry_view* view, size_t kmer_size,
                                         size_t start_pos, stack_type stack) noexcept
    : _view{ view }
      , _kmer_size{ kmer_size }
      , _start_pos{ start_pos }
      , _stack{ std::move(stack) }
      , _threshold{ powf(1.0f / (view->get_end_pos() - view->get_start_pos()), view->get_entry()->get_alphabet_size()) }
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
    calculate_next_phylokmer();
    return *this;
}

phylo_kmer_iterator::reference phylo_kmer_iterator::operator*()
{
    return _stack.back().mmer;
}

phylo_kmer_iterator::pointer phylo_kmer_iterator::operator->()
{
    return &(_stack.back().mmer);
}

void phylo_kmer_iterator::calculate_next_phylokmer()
{
    auto has_ended = bottom_up();
    if (has_ended)
    {
        return;
    }
    top_down();
}

bool phylo_kmer_iterator::bottom_up()
{
    const auto alphabet_size = _view->get_entry()->get_alphabet_size();

    /// go bottom-up and pop out m-mers (m <= k)
    while (!_stack.empty() && _stack.back().last_index >= alphabet_size - 1)
    {
        _stack.pop_back();
    }

    /// if we discovered the last k-mer, iterator has ended
    if (_stack.empty())
    {
        *this = _view->end();
        return true;
    }
    return false;
}

void phylo_kmer_iterator::top_down()
{
    const auto alphabet_size = _view->get_entry()->get_alphabet_size();
    const auto kmask = mask(_kmer_size, alphabet_size);

    auto last_mmer = _stack.back();
    auto next_row = last_mmer.last_row;
    auto next_index = last_mmer.last_index + 1;
    const auto& last_letter = _view->at(last_mmer.last_row, last_mmer.last_index);
    auto next_kmer_value = right_shift(last_mmer.mmer.value, _kmer_size, alphabet_size);
    auto next_kmer_score = last_mmer.mmer.score - last_letter.score;
    _stack.pop_back();

    /// go top-down and push undiscovered m-mers (m <= k)
    while (_stack.size() < _kmer_size)
    {
        /// calculate new k-mer value and score
        const auto& new_letter = _view->at(next_row, next_index);
        next_kmer_value = (left_shift(next_kmer_value, _kmer_size, alphabet_size) & kmask) | new_letter.index;
        next_kmer_score += new_letter.score;
        _stack.push_back(phylo_mmer{ phylo_kmer{ next_kmer_value, next_kmer_score }, next_row, next_index });

        ++next_row;
        next_index = 0;
    }
}


branch_entry_view::branch_entry_view(const branch_entry* entry, size_t start, size_t end) noexcept
    : _entry{ entry }
      , _start{ start }
      , _end{ end }
{}

branch_entry_view::branch_entry_view(const branch_entry_view& other) noexcept
    : branch_entry_view{ other._entry, other._start, other._end}
{}

branch_entry_view::const_iterator branch_entry_view::begin() const
{
    return phylo_kmer_iterator{ this, _end - _start, _start, populate_stack() };
}

branch_entry_view::const_iterator branch_entry_view::end() const
{
    return phylo_kmer_iterator{ this, 0, 0, {}};
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

phylo_kmer_iterator::stack_type branch_entry_view::populate_stack() const
{
    const auto kmer_size = _end - _start;
    const auto alphabet_size = _entry->get_alphabet_size();

    auto stack = phylo_kmer_iterator::stack_type{};
    phylo_kmer kmer { 0, 0.0f };
    for (size_t i = 0; i < kmer_size; ++i)
    {
        const size_t next_index = 0;
        const auto& proba_pair = _entry->at(i, next_index);
        kmer.value = left_shift(kmer.value, kmer_size, alphabet_size) | proba_pair.index;
        kmer.score += proba_pair.score;
        stack.push_back(phylo_kmer_iterator::phylo_mmer{ kmer, i, next_index });
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
