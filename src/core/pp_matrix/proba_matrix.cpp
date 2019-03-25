#include "proba_matrix.h"
#include "utils/algorithm.h"
#include <algorithm>

size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return std::begin(_data)->second.alignment_size();
}

size_t proba_matrix::num_variants() const
{
    return std::begin(_data)->second.alphabet_size();
}

proba_matrix::mapped_type& proba_matrix::operator[](key_t branch_id)
{
    return _data[branch_id];
}

const proba_matrix::mapped_type& proba_matrix::at(key_t branch_id) const
{
    return _data.at(branch_id);
}

proba_matrix::iterator proba_matrix::find(const key_t& key)
{
    return _data.find(key);
}

proba_matrix::const_iterator proba_matrix::find(const key_t& key) const
{
    return _data.find(key);
}

proba_matrix::iterator proba_matrix::begin()
{
    return std::begin(_data);
}

proba_matrix::iterator proba_matrix::end()
{
    return std::end(_data);
}

proba_matrix::const_iterator proba_matrix::begin() const
{
    return _data.cbegin();
}

proba_matrix::const_iterator proba_matrix::end() const
{
    return std::end(_data);
}
