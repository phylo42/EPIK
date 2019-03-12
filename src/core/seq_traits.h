#pragma once

#include <set>

class seq_traits
{
public:
    seq_traits(const std::set<char>& char_set, const std::set<char>& ambiguous_chars);
    seq_traits(const seq_traits&) = default;
    seq_traits(seq_traits&&) = default;
    ~seq_traits() = default;

    bool is_valid(char c) const;
    bool is_ambiguous(char c) const;
    bool is_gap(char c) const;
private:
    const std::set<char> _char_set;
    const std::set<char> _ambiguous_chars;
};

seq_traits make_dna_seq_traits();