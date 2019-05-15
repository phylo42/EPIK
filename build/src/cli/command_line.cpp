#include <set>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "command_line.h"
#include "../return.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

void check_conflicts(const po::variables_map& vm,
                     const std::string& opt1, const std::string& opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() &&
        vm.count(opt2) && !vm[opt2].defaulted())
    {
        //throw conflicting_options(opt1, opt2);
        throw std::runtime_error("Conflicting options: " + opt1 + ", " + opt2);
    }
}

//--------------------------------------------------------------------------
namespace cli
{
    static std::string HELP = "help", HELP_SHORT = "h";
    static std::string WORKING_DIR = "workdir", WORKING_DIR_SHORT = "w";
    static std::string AR_PROBABILITIES = "ar-probabilities", AR_PROBABILITIES_SHORT = "a";
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";
    static std::string EXTENDED_TREE = "extended_tree", EXTENDED_TREE_SHORT = "x";
    static std::string EXTENDED_MAPPING = "extended_mapping", EXTENDED_MAPPING_SHORT = "e";
    static std::string ARTREE_MAPPING = "artree_mapping", ARTREE_MAPPING_SHORT = "m";
    static std::string K = "k", K_SHORT = "k";

    const po::options_description get_opt_description()
    {
        po::options_description desc("General options");
        desc.add_options()
            ((HELP + "," + HELP_SHORT).c_str(),
             "Show help")
            ((WORKING_DIR + "," + WORKING_DIR_SHORT).c_str(), po::value<fs::path>()->default_value(fs::current_path()),
             "Path to the working directory")
            ((AR_PROBABILITIES + "," + AR_PROBABILITIES_SHORT).c_str(), po::value<fs::path>()->required(),
             "Ancestral reconstruction probabilities file")
            ((REFTREE + "," + REFTREE_SHORT).c_str(), po::value<fs::path>()->required(),
             "Original phylogenetic tree file")
            ((EXTENDED_TREE + "," + EXTENDED_TREE_SHORT).c_str(), po::value<fs::path>()->required(),
             "Extended phylogenetic tree file")
            ((EXTENDED_MAPPING + "," + EXTENDED_MAPPING_SHORT).c_str(), po::value<fs::path>()->required(),
             "Original mapping file")
            ((ARTREE_MAPPING + "," + ARTREE_MAPPING_SHORT).c_str(), po::value<fs::path>()->required(),
             "Ancestral reconstruction tree mapping file")
            ((K + "," + K_SHORT).c_str(), po::value<size_t>()->default_value(8),
             "k-mer length used at DB build");
        return desc;
    }

    const std::string get_option_list()
    {
        std::stringstream ss;
        ss << get_opt_description();
        return ss.str();
    }

    const cli_parameters process_command_line(int argc, const char* argv[])
    {
        cli_parameters parameters;
        try
        {
            const po::options_description desc = get_opt_description();
            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv, desc), vm);
            po::notify(vm);

            if (vm.count(HELP))
            {
                parameters.action = action_t::help;
                return parameters;
            }
            else
            {
                parameters.action = action_t::build;
            }

            parameters.working_directory = vm[WORKING_DIR].as<fs::path>().string();
            parameters.ar_probabilities_file = vm[AR_PROBABILITIES].as<fs::path>().string();
            parameters.original_tree_file = vm[REFTREE].as<fs::path>().string();
            parameters.extended_tree_file = vm[EXTENDED_TREE].as<fs::path>().string();
            parameters.extended_mapping_file = vm[EXTENDED_MAPPING].as<fs::path>().string();
            parameters.artree_mapping_file = vm[ARTREE_MAPPING].as<fs::path>().string();
            parameters.kmer_size = vm[K].as<size_t>();
        }
        catch (const po::error& e)
        {
            //throw bad_options(e.what());
            throw std::runtime_error(e.what());
        }
        return parameters;
    }
}

