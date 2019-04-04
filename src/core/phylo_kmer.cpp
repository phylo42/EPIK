#include "phylo_kmer.h"
#include <limits>
#include <type_traits>
#include <cmath>

/// Compares two floats for almost equality.
/// From: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp = 1) noexcept
{
    // The machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           // unless the result is subnormal
           || std::abs(x - y) < std::numeric_limits<T>::min();
}

const kmer_t nan_value = 0;
const score_t nan_score = std::numeric_limits<score_t>::quiet_NaN();

bool phylo_kmer::is_nan() const
{
    return (value == nan_value) && (score == nan_score);
}

bool operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept
{
    if (lhs.is_nan() || rhs.is_nan())
    {
        return false;
    }
    else
    {
        return (lhs.value == rhs.value) && (almost_equal<score_t>(lhs.score, rhs.score));
    }
}

phylo_kmer make_napk()
{
    return phylo_kmer { nan_value, nan_score };
}

score_t score_threshold(size_t kmer_size)
{
    return std::log10(powf(1.0f / seq_traits<seq_type>::alphabet_size, float(kmer_size)));
}

