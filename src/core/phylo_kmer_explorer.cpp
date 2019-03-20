#include "phylo_kmer_explorer.h"

phylo_kmer_explorer::phylo_kmer_explorer(const seq_traits& traits, const proba_matrix::branch_entry_t& probas)
    : _seq_traits(traits)
    , _probas(probas)
{}