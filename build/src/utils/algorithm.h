#ifndef RAPPAS_CPP_ALGORITHM_H
#define RAPPAS_CPP_ALGORITHM_H

#include <cstddef>
#include <iterator>
#include <vector>
#include <numeric>
#include <algorithm>

/// Auxiliary compare class for the argsort function
template<class RandomInputIt, class RandomOutputIt, class Compare>
struct index_compare
{
    /// WARNING:
    /// Actually, RandomOutputIt::difference_type and RandomInputIt::difference_type must be
    /// the same, because we have to use the difference type to address the data
    /// and the index structures at the same time
    using output_diff_t = typename std::iterator_traits<RandomOutputIt>::difference_type;

    explicit index_compare(RandomOutputIt idx_first, RandomInputIt data_first, Compare compare)
        : _idx_first(idx_first)
        , _data_first(data_first)
        , _compare(compare)
    {}

    bool operator()(std::ptrdiff_t i, std::ptrdiff_t j) const
    {
        /// WARNING: this expression can be ill-typed in the case when then RandomOutputIt::difference_type
        /// is not compatible with the RandomInputIt::difference_type. In practice, it is not the case if we use
        /// STL containers (which usually are implemented with std::ptrdiff_t as an iterator difference type.
        /// See more: https://en.cppreference.com/w/cpp/iterator/iterator_traits
        return _compare(*(_data_first + i), *(_data_first + j));
    }

private:
    RandomOutputIt _idx_first;
    RandomInputIt _data_first;
    Compare _compare;
};


/// Argsort function. Sorts the indexes of input array (the input array remains unsorted)
template<class RandomInputIt, class RandomOutputIt, class Compare = std::less<>>
void argsort(RandomOutputIt first, RandomOutputIt last,
             RandomInputIt data_first,
             Compare comp = Compare())
{
    std::sort(first, last, index_compare<RandomInputIt, RandomOutputIt, Compare>(first, data_first, comp));
}

#endif //RAPPAS_CPP_ALGORITHM_H
