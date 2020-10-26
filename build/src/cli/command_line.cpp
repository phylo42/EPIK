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
    static std::string REFALIGN = "refalign";
    static std::string AR_PROBABILITIES = "ar-probabilities", AR_PROBABILITIES_SHORT = "a";
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";

    static std::string AR_MODEL = "model";
    static std::string REDUCTION_RATIO = "reduction-ratio";
    static std::string EXTENDED_TREE = "extended-tree", EXTENDED_TREE_SHORT = "x";
    static std::string EXTENDED_MAPPING = "extended-mapping", EXTENDED_MAPPING_SHORT = "e";
    static std::string ARTREE_MAPPING = "artree-mapping", ARTREE_MAPPING_SHORT = "m";
    static std::string K = "k", K_SHORT = "k";
    static std::string OMEGA="omega", OMEGA_SHORT="o";
    static std::string NUM_THREADS = "num_threads", NUM_THREADS_SHORT = "j";

    // Filtering options
    static std::string MU = "mu", MU_SHORT = "u";
    static std::string NO_FILTER = "no-filter";
    static std::string ENTROPY = "entropy";
    static std::string MAX_DEVIATION = "max-deviation";
    static std::string LOG_MAX_DEVIATION = "log-max-deviation";
    static std::string MAX_DIFF = "max-difference";
    static std::string LOG_MAX_DIFF = "log-max-difference";
    static std::string STD_DEVIATION = "sd";
    static std::string LOG_STD_DEVIATION = "log-sd";
    static std::string RANDOM = "random";
    static std::string MERGE_BRANCHES = "merge-branches";

    bool no_filter_flag = true;
    bool entropy_flag = false;
    bool max_deviation_filter_flag = false;
    bool log_max_deviation_filter_flag = false;
    bool max_difference_filter_flag = false;
    bool log_max_difference_filter_flag = false;
    bool std_deviation_filter_flag = false;
    bool log_std_deviation_filter_flag = false;
    bool random_filter_flag = false;
    bool merge_branches_flag = false;

    po::options_description get_opt_description()
    {
        po::options_description desc("General options");
        desc.add_options()
            ((HELP + "," + HELP_SHORT).c_str(),
             "Show help")
            ((WORKING_DIR + "," + WORKING_DIR_SHORT).c_str(), po::value<fs::path>()->default_value(fs::current_path()),
             "Path to the working directory")
            (REFALIGN .c_str(), po::value<fs::path>()->required(),
             "Reference alignment in fasta format."
             "It must be the multiple alignment from which the reference tree was built.")
            ((AR_PROBABILITIES + "," + AR_PROBABILITIES_SHORT).c_str(), po::value<fs::path>()->required(),
             "Ancestral reconstruction probabilities file")
            ((REFTREE + "," + REFTREE_SHORT).c_str(), po::value<fs::path>()->required(),
             "Original phylogenetic tree file")
            (AR_MODEL.c_str(), po::value<std::string>()->required(),
             "Model used in AR, one of the following:"
             "nucl  : JC69, HKY85, K80, F81, TN93, GTR"
             "amino : LG, WAG, JTT, Dayhoff, DCMut, CpREV, mMtREV, MtMam, MtArt")
            (REDUCTION_RATIO.c_str(), po::value<double>()->default_value(0.99),
             "Ratio for alignment reduction, e.g. sites holding >X% gaps are ignored.")
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
            ((MERGE_BRANCHES).c_str(), po::bool_switch(&merge_branches_flag))
            ((NO_FILTER).c_str(), po::bool_switch(&no_filter_flag))
            ((ENTROPY).c_str(), po::bool_switch(&entropy_flag))
            ((MAX_DEVIATION).c_str(), po::bool_switch(&max_deviation_filter_flag))
            ((LOG_MAX_DEVIATION).c_str(), po::bool_switch(&log_max_deviation_filter_flag))
            ((MAX_DIFF).c_str(), po::bool_switch(&max_difference_filter_flag))
            ((LOG_MAX_DIFF).c_str(), po::bool_switch(&log_max_difference_filter_flag))
            ((STD_DEVIATION).c_str(), po::bool_switch(&std_deviation_filter_flag))
            ((LOG_STD_DEVIATION).c_str(), po::bool_switch(&log_std_deviation_filter_flag))
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
            parameters.alignment_file = vm[REFALIGN].as<fs::path>().string();
            parameters.ar_probabilities_file = vm[AR_PROBABILITIES].as<fs::path>().string();
            parameters.original_tree_file = vm[REFTREE].as<fs::path>().string();
            parameters.ar_model = vm[AR_MODEL].as<std::string>();
            parameters.reduction_ratio = vm[REDUCTION_RATIO].as<double>();
            parameters.extended_tree_file = vm[EXTENDED_TREE].as<fs::path>().string();
            parameters.extended_mapping_file = vm[EXTENDED_MAPPING].as<fs::path>().string();
            parameters.artree_mapping_file = vm[ARTREE_MAPPING].as<fs::path>().string();
            parameters.kmer_size = vm[K].as<size_t>();
            parameters.omega = vm[OMEGA].as<xpas::phylo_kmer::score_type>();
            parameters.num_threads = vm[NUM_THREADS].as<size_t>();
            parameters.mu = vm[MU].as<double>();
            parameters.no_filter = no_filter_flag;
            parameters.entropy_filter = entropy_flag;
            parameters.max_dev_filter = max_deviation_filter_flag;
            parameters.log_max_dev_filter = log_max_deviation_filter_flag;
            parameters.max_diff_filter = max_difference_filter_flag;
            parameters.log_max_diff_filter = log_max_difference_filter_flag;
            parameters.random_filter = random_filter_flag;
            parameters.std_dev_filter = std_deviation_filter_flag;
            parameters.log_std_dev_filter = log_std_deviation_filter_flag;
            parameters.merge_branches = merge_branches_flag;
        }
        catch (const po::error& e)
        {
            throw std::runtime_error(e.what());
        }
        return parameters;
    }
}

