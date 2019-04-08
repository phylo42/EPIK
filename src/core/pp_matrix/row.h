#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <core/phylo_kmer.h>
#include <array>

struct proba_pair
{
    score_t score;
    size_t index;
};

/// alignment position type
using pos_t = size_t;
constexpr pos_t not_a_position = std::numeric_limits<pos_t>::max();

using branch_id = uint16_t;
using row = std::array<proba_pair, seq_traits<seq_type>::alphabet_size>;

#endif //RAPPAS_CPP_ROW_H
