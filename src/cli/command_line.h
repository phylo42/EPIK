#pragma once

#include <exception>
#include <string>
#include <map>

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
        action_t action;
        std::string tree_file;
        std::string working_directory;

        // output options
        // ...

        // algo options
        size_t kmer_size;
        double reduction_ratio;

        // debug options
        std::string ar_directory;

        // TODO: support all options
    };

    const std::string get_option_list();
    const cli_parameters process_command_line(int argc, const char* argv[]);
}