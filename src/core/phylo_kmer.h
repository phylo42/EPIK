#pragma once

#include <cstdint>

/// TODO: Remove redundant fields from here. This structure is used in two different places: phylokmer calculation
/// and as an entry type for a phylokmer database. Split it
struct phylo_kmer
{
    using kmer_value_t = uint32_t;
    using branch_node_t = uint16_t;

    kmer_value_t kmer_value;
    branch_node_t branch_node;
    float score;
};