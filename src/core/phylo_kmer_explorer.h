#pragma once

#include <unordered_map>
#include "seq_traits.h"
#include "proba_matrix.h"
#include "phylo_kmer.h"

class phylo_kmer_iterator
{
public:
    // TODO:
    // Sorting of columns here
    // store sorted indexes in proba_matrix?
};

class phylo_kmer_explorer
{
public:
    phylo_kmer_explorer(const seq_traits& traits, const proba_matrix::branch_entry_t& probas);

private:
    seq_traits _seq_traits;
    const proba_matrix::branch_entry_t& _probas;
};

