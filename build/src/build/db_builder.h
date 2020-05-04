#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>

namespace xpas
{
    class phylo_kmer_db;
}

namespace rappas
{
    xpas::phylo_kmer_db build(const std::string& working_directory, const std::string& ar_probabilities_file,
        const std::string& original_tree_file, const std::string& extended_tree_file,
        const std::string& extended_mapping_file, const std::string& artree_mapping_file,
        size_t kmer_size, xpas::phylo_kmer::score_type omega, size_t num_threads);
}

#endif