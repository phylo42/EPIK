#include <algorithm>
#include "seq_traits.h"

using std::vector;
using char_t = seq_traits::char_type;

seq_traits::seq_traits(const vector<char_t>& char_set, const vector<char_t>& ambiguous_chars)
    : char_set{ char_set }
    , ambiguous_chars{ ambiguous_chars }
{}

bool seq_traits::is_valid(char_t c) const
{
    return find(begin(char_set), end(char_set), c) != end(char_set);
}

bool seq_traits::is_ambiguous(char_t c) const
{
    return find(begin(ambiguous_chars), end(ambiguous_chars), c) != end(ambiguous_chars);
}

bool seq_traits::is_gap(char_t c) const
{
    return c == '-' || c == '.';
}

seq_traits make_dna_seq_traits()
{
    const vector<char_t> char_set = { 'A', 'C', 'G', 'T' };
    const vector<char_t> ambiguous_chars = {
        'N', '.', '-',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
    };
    return seq_traits{ char_set, ambiguous_chars };
}