#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <core/phylo_kmer.h>
#include <array>

namespace rappas
{
    struct proba_pair
    {
        core::phylo_kmer::score_type score;
        size_t index;
    };

    using branch_type = core::phylo_kmer::branch_type;
    using row_type = std::array<proba_pair, core::seq_traits::alphabet_size>;
}

#endif //RAPPAS_CPP_ROW_H
