#ifndef RAPPAS_CPP_SEQ_TRAITS_H
#define RAPPAS_CPP_SEQ_TRAITS_H

#include <set>

struct seq_traits
{
    seq_traits(const std::set<char>& char_set, const std::set<char>& ambiguous_chars);
    seq_traits(const seq_traits&) = default;
    seq_traits(seq_traits&&) = default;
    seq_traits& operator=(const seq_traits&) = default;
    seq_traits& operator=(seq_traits&&) = default;
    ~seq_traits() noexcept = default;

    bool is_valid(char c) const;
    bool is_ambiguous(char c) const;
    bool is_gap(char c) const;

    std::set<char> char_set;
    std::set<char> ambiguous_chars;
};

static const auto dna_seq_traits = seq_traits {
    {'A', 'C', 'G', 'T'},
    {'N', '.', '-',
     'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}
};

#endif