#ifndef RAPPAS_CPP_SEQ_TRAITS_H
#define RAPPAS_CPP_SEQ_TRAITS_H

#include <vector>
#include <cstdint>

struct seq_traits
{
    using char_type = uint8_t;

    seq_traits(const std::vector<char_type>& char_set, const std::vector<char_type>& ambiguous_chars);
    seq_traits(const seq_traits&) = default;
    seq_traits(seq_traits&&) = default;
    seq_traits& operator=(const seq_traits&) = default;
    seq_traits& operator=(seq_traits&&) = default;
    ~seq_traits() noexcept = default;

    bool is_valid(char_type c) const;
    bool is_ambiguous(char_type c) const;
    bool is_gap(char_type c) const;

    std::vector<char_type> char_set;
    std::vector<char_type> ambiguous_chars;
};

static const auto dna_seq_traits = seq_traits {
    {'A', 'C', 'G', 'T'},
    {'N', '.', '-',
     'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'}
};

#endif