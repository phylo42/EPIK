#include "seq_traits.h"

seq_traits::seq_traits(const std::set<char>& char_set, const std::set<char>& ambiguous_chars)
    : char_set(char_set)
     , ambiguous_chars(ambiguous_chars)
{}

bool seq_traits::is_valid(char c) const
{
    return char_set.find(c) != char_set.end();
}

bool seq_traits::is_ambiguous(char c) const
{
    return ambiguous_chars.find(c) != ambiguous_chars.end();
}

bool seq_traits::is_gap(char c) const
{
    return c == '-' || c == '.';
}

seq_traits make_dna_seq_traits()
{
    const std::set<char> char_set = {'A', 'C', 'G', 'T'};
    const std::set<char> ambiguous_chars = {
        'N', '.', '-',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
    };
    return seq_traits(char_set, ambiguous_chars);
}