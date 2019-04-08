#ifndef RAPPAS_CPP_SEQ_H
#define RAPPAS_CPP_SEQ_H

#include <cstdint>
#include <cstddef>

struct dna
{};

typedef uint64_t kmer_t;

using seq_type = dna;

template<typename SeqType>
struct seq_traits;

template<>
struct seq_traits<dna>
{
    using char_type = uint8_t;
    static constexpr char_type char_set[] = {'A', 'C', 'G', 'T'};
    static constexpr char_type ambiguous_chars[] =  {'N', '.', '-', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
    static constexpr size_t alphabet_size = sizeof(char_set);
    static constexpr size_t max_kmer_length = 16;
};



template<typename SeqType>
constexpr kmer_t bit_length();

template<>
constexpr kmer_t bit_length<dna>()
{
    return 2ul;
}

template<typename SeqType>
constexpr kmer_t rightest_symbol_mask();

template<>
constexpr kmer_t rightest_symbol_mask<dna>()
{
    return ~0b11ul;
}

#endif