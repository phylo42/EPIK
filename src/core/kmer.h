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
/// size alphabet_size.
/// TODO: implement in non-recursively
inline kmer_t mask(size_t kmer_size, size_t alphabet_size)
{
    if (kmer_size == 1)
    {
        return 1;
    }
    return (mask(kmer_size - 1, alphabet_size) << bits_required(1, alphabet_size)) | mask(1, alphabet_size);
};

inline kmer_t shift_append(kmer_t value, size_t kmer_size, size_t alphabet_size, unsigned char letter)
{
    return ((value << bits_required(1, alphabet_size)) | letter);
}

#endif