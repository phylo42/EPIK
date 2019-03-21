#include <algorithm>
#include "proba_matrix.h"
#include "../utils/algorithm.h"

proba_matrix::proba_matrix()
    : _data()
{}

proba_matrix::proba_matrix(proba_matrix&& other)
    : _data(std::move(other._data))
{}

size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return begin(_data)->second.size();
}

size_t proba_matrix::num_variants() const
{
    return begin(_data)->second.begin()->first.size();
}

const proba_matrix::branch_entry_t proba_matrix::at(int branch_id) const
{
    return _data.at(branch_id);
}

void proba_matrix::add_branch_entry(int node_id, const branch_entry_t& branch_entry)
{
    _data.emplace(node_id, branch_entry);
}

void proba_matrix::sort()
{
    for (auto& map_entry : _data)
    {
        for (auto& row : map_entry.second)
        {
            argsort(begin(row.second), end(row.second), begin(row.first), std::greater<>());
        }
    }
}
