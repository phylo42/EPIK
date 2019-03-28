#ifndef RAPPAS_CPP_PHYLO_KMER_H
#define RAPPAS_CPP_PHYLO_KMER_H

#include <cstdint>

#include "kmer.h"

typedef float score_t;
typedef uint16_t branch_node_t;

struct phylo_db_entry
{
    kmer_t value;
    score_t score;
    branch_node_t branch_node;
};

struct phylo_kmer
{
    bool is_nan() const;

    kmer_t value;
    score_t score;
};

bool operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept;

/// Returns a phylo_kmer with special values, considered as NotAPhyloKmer. This phylo kmer
/// can not be equeal to any other phylo kmer (including itself)
phylo_kmer make_napk();

#endif