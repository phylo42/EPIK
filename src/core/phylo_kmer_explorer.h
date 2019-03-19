#pragma once

#include <unordered_map>
#include "seq_traits.h"
#include "ar.h"
#include "phylo_kmer.h"


class phylo_kmer_explorer
{
public:
    phylo_kmer_explorer(const seq_traits& traits, proba_matrix&& probas);

    bool has_next() const;

    phylo_kmer next_phylo_kmer();
private:
    seq_traits _seq_traits;
    proba_matrix _probas;
};

