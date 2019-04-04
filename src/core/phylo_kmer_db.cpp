#include "phylo_kmer_db.h"

void phylo_kmer_db::put(kmer_t key, branch_node_t branch, score_t score, size_t position)
{
    if (auto it = _map.find(key); it != _map.end())
    {
        if (auto inner_map_it = it->second.find(branch); inner_map_it != it->second.end())
        {
            if (inner_map_it->second.score < score)
            {
                _map[key][branch] = { score, position };
            }
        }
        else
        {
            _map[key][branch] = { score, position };
        }
    }
    else
    {
        _map[key][branch] = { score, position };
    }
}


phylo_kmer_db::const_iterator phylo_kmer_db::begin() const
{
    return std::begin(_map);
}

phylo_kmer_db::const_iterator phylo_kmer_db::end() const
{
    return std::end(_map);
}

size_t phylo_kmer_db::size() const
{
    return _map.size();
}

