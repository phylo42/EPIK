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
        std::string ar_probabilities_file;
        std::string original_tree_file;
        std::string extended_tree_file;
        std::string extended_mapping_file;
        std::string artree_mapping_file;

        // algo options
        size_t kmer_size;
        xpas::phylo_kmer::score_type omega;
        size_t num_threads;

        // k-mer filtering paramaters, mutually exclusive
        bool no_filter;
        bool entropy_filter;
        bool maxdev_filter;
        bool maxdiff_filter;
        bool random_filter;
        bool std_dev_filter;

        // k-mer filtering threshold
        double mu;
    };

    std::string get_option_list();
    cli_parameters process_command_line(int argc, const char* argv[]);
}

#endif