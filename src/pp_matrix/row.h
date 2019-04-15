#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <core/phylo_kmer.h>
#include <array>

struct alignment
{
    /// alignment position type
    using pos_type = size_t;

    static constexpr pos_type not_a_position = std::numeric_limits<pos_type>::max();
};

struct proba_pair
{
    core::phylo_kmer::score_type score;
    size_t index;
};

using branch_id = uint16_t;
using row = std::array<proba_pair, core::seq_traits::alphabet_size>;

#endif //RAPPAS_CPP_ROW_H
