#include <iostream>
#include <fstream>
#include "cli/command_line.h"
#include "cli/exceptions.h"
#include "core/return.h"
#include "core/db_builder.h"

return_code_t print_help()
{
    std::cout << "RAPPAS2" << std::endl << std::endl
              << "Usage: rappas2 [...]" << std::endl
              << cli::get_option_list() << std::endl;

    return return_code::help;
}

return_code_t run(const cli::cli_parameters& parameters)
{
    switch (parameters.action)
    {
        case cli::help:
        {
            return print_help();
        }
        case cli::build:
        {
            db_builder builder(parameters.working_directory, parameters.ar_probabilities_file,
                               parameters.tree_file, parameters.extended_mapping_file, parameters.artree_mapping_file,
                               parameters.kmer_size, dna_seq_traits);
            return builder.run();
        }
        default:
        {
            return return_code::unknown_error;
        }
    };
}

int main(int argc, const char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        const cli::cli_parameters parameters = cli::process_command_line(argc, argv);
        return run(parameters);
    }
    catch (const conflicting_options& e)
    {
        std::cout << e.what() << std::endl;
        return return_code::argument_error;
    }
    catch (const bad_options& e)
    {
        std::cout << e.what() << std::endl;
        return return_code::argument_error;
    }
    catch (const std::runtime_error& e)
    {
        std::cout << e.what() << std::endl;
        return return_code::argument_error;
    }
    /*catch (...)
    {
        std::cout << "Unexpected error. " << std::endl;
        return return_code::unknown_error;
    }*/
}