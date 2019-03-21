#include "proba_matrix.h"
#include "../utils/algorithm.h"

#include <algorithm>


size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return std::begin(_data)->second.size();
}

size_t proba_matrix::num_variants() const
{
    return std::begin(_data)->second.begin()->first.size();
}

proba_matrix::branch_entry_t& proba_matrix::operator[](key_t branch_id)
{
    return _data[branch_id];
}

const proba_matrix::branch_entry_t proba_matrix::at(key_t branch_id) const
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

void proba_matrix::sort()
{
    for (auto& map_entry : _data)
    {
        for (auto& row : map_entry.second)
        {
            argsort(std::begin(row.second), std::end(row.second),
                    std::begin(row.first),
                    std::greater<>());
        }
    }
}
