#pragma once

#include <exception>
#include <string>
#include <map>
#include "../ar/evo_model.h"

namespace cli
{
    enum action_t
    {
        build = 0,
        place = 1,
        help = 2
    };

    struct cli_parameters
    {
        // main options
        std::string ar_binary;
        action_t action;
        std::string alignment_file;
        std::string tree_file;
        std::string working_directory;

        // output options
        // ...

        // algo options
        ar::evo_model model = ar::get_default_model();
        std::string ar_parameters;
        size_t kmer_size;
        double reduction_ratio;

        // debug options
        std::string ar_directory;

        // TODO: support all options
    };

    const std::string get_option_list();
    const cli_parameters process_command_line(int argc, const char* argv[]);
}