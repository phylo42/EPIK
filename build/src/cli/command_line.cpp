#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "command_line.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

//--------------------------------------------------------------------------
namespace cli
{
    static std::string HELP = "help", HELP_SHORT = "h";
    static std::string WORKING_DIR = "workdir", WORKING_DIR_SHORT = "w";
    static std::string AR_PROBABILITIES = "ar-probabilities", AR_PROBABILITIES_SHORT = "a";
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";
    static std::string EXTENDED_TREE = "extended-tree", EXTENDED_TREE_SHORT = "x";
    static std::string EXTENDED_MAPPING = "extended-mapping", EXTENDED_MAPPING_SHORT = "e";
    static std::string ARTREE_MAPPING = "artree-mapping", ARTREE_MAPPING_SHORT = "m";
    static std::string K = "k", K_SHORT = "k";
    static std::string OMEGA="omega", OMEGA_SHORT="o";
    static std::string NUM_THREADS = "num_threads", NUM_THREADS_SHORT = "j";
    static std::string MU = "mu", MU_SHORT = "u";
    static std::string NO_FILTER = "no-filter";
    static std::string ENTROPY = "entropy";
    static std::string MAX_DEVIATION = "max-deviation";
    static std::string MAX_DIFF = "max-difference";
    static std::string STD_DEVIATION = "sd";
    static std::string RANDOM = "random";

    bool no_filter_flag = true;
    bool entropy_flag = false;
    bool max_deviation_filter_flag = false;
    bool max_difference_filter_flag = false;
    bool random_filter_flag = false;
    bool std_deviation_filter_flag = false;

    po::options_description get_opt_description()
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
             "k-mer length used at DB build")
            ((OMEGA + "," + OMEGA_SHORT).c_str(), po::value<xpas::phylo_kmer::score_type>()->default_value(1.5),
             "Score threshold parameter")
            ((NUM_THREADS + "," + NUM_THREADS_SHORT).c_str(), po::value<size_t>()->default_value(1),
             "Number of threads")
            ((NO_FILTER).c_str(), po::bool_switch(&no_filter_flag))
            ((ENTROPY).c_str(), po::bool_switch(&entropy_flag))
            ((MAX_DEVIATION).c_str(), po::bool_switch(&max_deviation_filter_flag))
            ((MAX_DIFF).c_str(), po::bool_switch(&max_difference_filter_flag))
            ((STD_DEVIATION).c_str(), po::bool_switch(&std_deviation_filter_flag))
            ((RANDOM).c_str(), po::bool_switch(&random_filter_flag))
            ((MU + "," + MU_SHORT).c_str(), po::value<double>()->default_value(0.8));
        return desc;
    }

    std::string get_option_list()
    {
        std::stringstream ss;
        ss << get_opt_description();
        return ss.str();
    }

    cli_parameters process_command_line(int argc, const char* argv[])
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
            parameters.omega = vm[OMEGA].as<xpas::phylo_kmer::score_type>();
            parameters.num_threads = vm[NUM_THREADS].as<size_t>();
            parameters.mu = vm[MU].as<double>();
            parameters.no_filter = no_filter_flag;
            parameters.entropy_filter = entropy_flag;
            parameters.maxdev_filter = max_deviation_filter_flag;
            parameters.maxdiff_filter = max_difference_filter_flag;
            parameters.random_filter = random_filter_flag;
            parameters.std_dev_filter = std_deviation_filter_flag;
        }
        catch (const po::error& e)
        {
            throw std::runtime_error(e.what());
        }
        return parameters;
    }
}

