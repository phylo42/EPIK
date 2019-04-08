#include "phylo_kmer_db.h"

score_pos::score_pos() noexcept
    : score_pos(0.0f, not_a_position)
{}

score_pos::score_pos(score_t s, pos_t p) noexcept
{}

score_pos::operator bool() const
{
    return position != not_a_position;
}

void phylo_kmer_db::put(kmer_t key, branch_node_t branch, score_t score, size_t position)
{
    if (auto it = _map.find(key); it != _map.end())
    {
        auto& branch_entry = it->second[branch];
        if (branch_entry)
        {
            if (branch_entry.score < score)
            {
                branch_entry = { score, position };
            }
        }
        else
        {
            branch_entry = { score, position };
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

