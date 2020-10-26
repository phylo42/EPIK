#ifndef RAPPAS_CPP_COMMAND_LINE_H
#define RAPPAS_CPP_COMMAND_LINE_H

#include <exception>
#include <string>
#include <map>
#include <xpas/phylo_kmer.h>

namespace cli
{
    enum action_t
    {
        build = 0,
        help = 2
    };

    struct cli_parameters
    {
        // main options
        action_t action;
        std::string working_directory;
        std::string alignment_file;
        std::string ar_probabilities_file;
        std::string original_tree_file;
        std::string extended_tree_file;
        std::string extended_mapping_file;
        std::string artree_mapping_file;

        std::string ar_model;

        // Alignment filtering options
        double reduction_ratio;
        bool no_reduction;

        // algo options
        size_t kmer_size;
        xpas::phylo_kmer::score_type omega;
        size_t num_threads;

        bool merge_branches;

        // k-mer filtering paramaters, mutually exclusive
        bool no_filter;
        bool entropy_filter;
        bool max_dev_filter;
        bool log_max_dev_filter;
        bool max_diff_filter;
        bool log_max_diff_filter;
        bool random_filter;
        bool std_dev_filter;
        bool log_std_dev_filter;


        // k-mer filtering threshold
        double mu;
    };

    std::string get_option_list();
    cli_parameters process_command_line(int argc, const char* argv[]);
}

#endif