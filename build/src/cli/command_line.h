#ifndef RAPPAS_CPP_COMMAND_LINE_H
#define RAPPAS_CPP_COMMAND_LINE_H

#include <exception>
#include <string>
#include <map>
#include <core/phylo_kmer.h>

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
        core::phylo_kmer::score_type omega;
        size_t num_threads;

        double mu;
    };

    const std::string get_option_list();
    const cli_parameters process_command_line(int argc, const char* argv[]);
}

#endif