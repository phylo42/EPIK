#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>

namespace xpas
{
    class phylo_kmer_db;
}

namespace rappas
{
    enum class filter_type
    {
        no_filter,
        entropy,
        max_deviation,
        log_max_deviation,
        max_difference,
        log_max_difference,
        standard_deviation,
        log_standard_deviation,
        random
    };

    xpas::phylo_kmer_db build(const std::string& working_directory, const std::string& ar_probabilities_file,
        const std::string& original_tree_file, const std::string& extended_tree_file,
        const std::string& extended_mapping_file, const std::string& artree_mapping_file,
        size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter, double mu, size_t num_threads);
}

#endif