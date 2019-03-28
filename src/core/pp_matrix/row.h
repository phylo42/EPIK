#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <vector>
#include <core/phylo_kmer.h>

struct proba_pair
{
    score_t score;
    size_t index;
};

using branch_id = uint16_t;
using row = std::vector<proba_pair>;

#endif //RAPPAS_CPP_ROW_H
