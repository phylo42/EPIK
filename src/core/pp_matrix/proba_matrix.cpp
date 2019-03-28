#include "proba_matrix.h"
#include "utils/algorithm.h"
#include <algorithm>

size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return std::begin(_data)->second.get_alignment_size();
}

size_t proba_matrix::num_variants() const
{
    return std::begin(_data)->second.get_alphabet_size();
}

proba_matrix::mapped_type& proba_matrix::operator[](branch_id id)
{
    return _data[id];
}

const proba_matrix::mapped_type& proba_matrix::at(branch_id id) const
{
    return _data.at(id);
}

proba_matrix::iterator proba_matrix::find(const branch_id& id)
{
    return _data.find(id);
}

proba_matrix::const_iterator proba_matrix::find(const branch_id& id) const
{
    return _data.find(id);
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
