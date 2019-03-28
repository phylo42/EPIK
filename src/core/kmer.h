#ifndef RAPPAS_CPP_KMER_H
#define RAPPAS_CPP_KMER_H

#include <cstddef>
#include <cmath>

typedef uint64_t kmer_t;

/// Returns how many bits required to write down a maximum numerical value of kmer_size chars
/// over alphabet of size alphabet_size
inline size_t bits_required(size_t kmer_size, size_t alphabet_size)
{
    return kmer_size * (size_t)ceil(alphabet_size / 2.0f);
};

/// Returns the maximum value that can be writen down with kmer_size chars over alphabet of
/// size alphabet_size
inline kmer_t mask(size_t kmer_size, size_t alphabet_size)
{
    return (1u << (kmer_size * bits_required(1, alphabet_size))) - 1;
};

inline kmer_t left_shift(kmer_t value, size_t kmer_size, size_t alphabet_size)
{
    return value << bits_required(1, alphabet_size);
}

inline kmer_t right_shift(kmer_t value, size_t kmer_size, size_t alphabet_size)
{
    return value >> bits_required(1, alphabet_size);
}


#endif