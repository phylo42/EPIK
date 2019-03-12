#include "seq_traits.h"

seq_traits::seq_traits(const std::set<char>& char_set, const std::set<char>& ambiguous_chars)
    : _char_set(char_set)
    , _ambiguous_chars(ambiguous_chars)
{}

bool seq_traits::is_valid(char c) const
{
    return _char_set.find(c) != _char_set.end();
}

bool seq_traits::is_ambiguous(char c) const
{
    return _ambiguous_chars.find(c) != _ambiguous_chars.end();
}

bool seq_traits::is_gap(char c) const
{
    return c == '-' || c == '.';
}

size_t seq_traits::charset_size() const
{
    return _char_set.size();
}


seq_traits make_dna_seq_traits()
{
    const std::set<char> char_set = { 'A', 'C', 'G', 'T', 'N', '.', '-' };
    const std::set<char> ambiguous_chars = {
            'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', '-',
    };
    return seq_traits(char_set, ambiguous_chars);
}