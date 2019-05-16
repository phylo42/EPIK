#include <iostream>
#include <boost/filesystem.hpp>
#include <chrono>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>
#include "cli/command_line.h"
#include "cli/exceptions.h"
#include "return.h"
#include "build/db_builder.h"

namespace fs = boost::filesystem;

return_code print_help()
{
    std::cout << "RAPPAS2" << std::endl << std::endl
              << "Usage: rappas2 [...]" << std::endl
              << cli::get_option_list() << std::endl;

    return return_code::help;
}

return_code run(const cli::cli_parameters& parameters)
{
    switch (parameters.action)
    {
        case cli::help:
        {
            return print_help();
        }
        case cli::build:
        {
            if (parameters.kmer_size > core::seq_traits::max_kmer_length)
            {
                std::cerr << "Maximum k-mer size allowed: " << core::seq_traits::max_kmer_length << std::endl;
                return return_code::argument_error;
            }

            const auto db = rappas::build(parameters.working_directory, parameters.ar_probabilities_file,
                                          parameters.original_tree_file, parameters.extended_tree_file,
                                          parameters.extended_mapping_file, parameters.artree_mapping_file,
                                          parameters.kmer_size);

            const auto db_filename = fs::path(parameters.working_directory) / "DB.union";

            std::cout << "Saving database to: " << db_filename.string() << "..." << std::endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            core::save(db, db_filename.string());
            std::cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << std::endl << std::endl;

            return return_code::success;
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
        run(parameters);
    }
    catch (const conflicting_options& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const bad_options& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unexpected error. " << std::endl;
        return 1;
    }
    return 0;
}