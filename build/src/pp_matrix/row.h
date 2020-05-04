#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <xpas/phylo_kmer.h>
#include <array>

namespace rappas
{
    struct proba_pair
    {
        xpas::phylo_kmer::score_type score;
        xpas::phylo_kmer::key_type index;
    };

    using branch_type = xpas::phylo_kmer::branch_type;
    using row_type = std::array<proba_pair, xpas::seq_traits::alphabet_size>;
}

#endif //RAPPAS_CPP_ROW_H
