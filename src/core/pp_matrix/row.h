#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <vector>
#include <core/phylo_kmer.h>

struct proba_pair
{
    proba_pair(score_t v, size_t i) : value { v }, index { i } { }
    proba_pair(const proba_pair&) = default;
    proba_pair& operator=(const proba_pair&) = default;
    ~proba_pair() noexcept = default;

    score_t value;
    size_t index;
};

using row = std::vector<proba_pair>;

#endif //RAPPAS_CPP_ROW_H
