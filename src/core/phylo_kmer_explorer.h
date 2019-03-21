#ifndef RAPPAS_CPP_PHYLO_KMER_EXPLORER_H
#define RAPPAS_CPP_PHYLO_KMER_EXPLORER_H

#include <unordered_map>
#include "seq_traits.h"
#include "proba_matrix.h"
#include "phylo_kmer.h"

class phylo_kmer_iterator
{
public:
};

template<size_t kmer_size>
class phylo_kmer_explorer
{
public:
    phylo_kmer_explorer(const seq_traits& traits, const proba_matrix::branch_entry_t& probas);
/*
    phylo_kmer_iterator begin()
    {
        auto branch_it = _probas.begin();
        uint64_t kmer_value = 0;
        for (size_t i = 0; i < kmer_size; ++i)
        {
            branch_it->first[branch_it->second];
        }
    }

    phylo_kmer_iterator end()
    {

    }*/

private:
    seq_traits _seq_traits;
    const proba_matrix::branch_entry_t& _probas;
};

#endif