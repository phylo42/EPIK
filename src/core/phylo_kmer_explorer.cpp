#include "phylo_kmer_explorer.h"



phylo_kmer_explorer::phylo_kmer_explorer(const seq_traits& traits, proba_matrix&& probas)
    : _seq_traits(traits)
    , _probas(std::move(probas))
{}


