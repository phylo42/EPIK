#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <xpas/phylo_kmer_db.h>

namespace rappas
{
    class alignment;

    enum class filter_type
    {
        no_filter,
        entropy,
        random
    };

    xpas::phylo_kmer_db build(std::string working_directory,
                              std::string ar_probabilities_file,
                              std::string original_tree_file,
                              std::string extended_tree_file,
                              std::string extended_mapping_file,
                              std::string artree_mapping_file,
                              rappas::alignment alignment,
                              bool merge_branches,
                              size_t kmer_size,
                              xpas::phylo_kmer::score_type omega,
                              filter_type filter,
                              double mu,
                              size_t num_threads);
}

#endif