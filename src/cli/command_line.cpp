#include <set>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "command_line.h"
#include "../core/return.h"

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
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";
    static std::string WORKING_DIR = "workdir", WORKING_DIR_SHORT = "w";
    static std::string K = "k", K_SHORT = "k";
    static std::string REDUCTION_RATIO = "ratio-reduction";

    void validate_reduction_ratio(double value)
    {
        if (value < 0.0 || value > 1.0)
        {
            throw po::validation_error(po::validation_error::invalid_option_value,
                    REDUCTION_RATIO, std::to_string(value));
        }
    }

    const po::options_description get_opt_description()
    {
        po::options_description desc("General options");
        desc.add_options()
                ((HELP + "," + HELP_SHORT).c_str(),
                        "Show help")
                ((REFTREE + "," + REFTREE_SHORT).c_str(), po::value<fs::path>()->required(),
                        "Phylogenetic tree file")
                ((WORKING_DIR + "," + WORKING_DIR_SHORT).c_str(), po::value<fs::path>()->default_value(fs::current_path()),
                        "Path to the working directory (b|p phase).")
                ((K + "," + K_SHORT).c_str(), po::value<size_t>()->default_value(8),
                         "k-mer length used at DB build. (b phase)")
                (REDUCTION_RATIO.c_str(), po::value<double>()->default_value(0.99f)->notifier(&validate_reduction_ratio),
                        "Ratio for alignment reduction, e.g. sites holding > X% gaps are ignored. (b phase)");
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

            if (vm.count(REFTREE))
            {
                parameters.tree_file = vm[REFTREE].as<fs::path>().string();
            }

            parameters.working_directory = vm[WORKING_DIR].as<fs::path>().string();
            parameters.kmer_size = vm[K].as<size_t>();
            parameters.reduction_ratio = vm[REDUCTION_RATIO].as<double>();
        }
        catch (const po::error& e)
        {
            //throw bad_options(e.what());
            throw std::runtime_error(e.what());
        }
        return parameters;
    }
}

