#ifndef RAPPAS_CPP_PHYLO_KMER_DB_H
#define RAPPAS_CPP_PHYLO_KMER_DB_H

#include <absl/container/flat_hash_map.h>
#include "phylo_kmer.h"

struct score_pos
{
    score_t score;
    size_t position;
};

class phylo_kmer_db
{
public:
    using storage = absl::flat_hash_map<kmer_t, absl::flat_hash_map<branch_node_t, score_pos>>;
    using const_iterator = storage::const_iterator;

    phylo_kmer_db() = default;
    phylo_kmer_db(const phylo_kmer_db&) = delete;
    phylo_kmer_db(phylo_kmer_db&&) = delete;
    phylo_kmer_db& operator=(const phylo_kmer_db&) = delete;
    phylo_kmer_db& operator=(phylo_kmer_db&&) = delete;
    ~phylo_kmer_db() noexcept = default;

    void put(kmer_t key, branch_node_t branch, score_t score, size_t position);

    const_iterator begin() const;
    const_iterator end() const;

    size_t size() const;
private:
    storage _map;
};


#endif
