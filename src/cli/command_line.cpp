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
    static std::string PHASE = "phase", PHASE_SHORT = "p";
    static std::string REFALIGN = "refalign", REFALIGN_SHORT = "r";
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";
    static std::string WORKING_DIR = "workdir", WORKING_DIR_SHORT = "w";
    static std::string ARDIR = "ardir", ARDIR_SHORT = "a";
    static std::string ARBINARY = "arbinary", ARBINARY_SHORT = "b";
    static std::string K = "k", K_SHORT = "k";
    static std::string MODEL = "model", MODEL_SHORT = "m";
    static std::string REDUCTION_RATIO = "ratio-reduction";

    struct phase_t
    {
        static std::string PHASE_BUILD, PHASE_PLACE;
        phase_t(const std::string& v) : value(v) { }
        std::string value;
    };
    std::string phase_t::PHASE_BUILD = "b";
    std::string phase_t::PHASE_PLACE = "p";

    void validate(boost::any& v, std::vector<std::string> const& values, phase_t*, int)
    {
        // Make sure no previous assignment to 'v' was made.
        po::validators::check_first_occurrence(v);

        // Extract the first string from 'values'. If there is more than
        // one string, it's an error, and exception will be thrown.
        std::string const& s = po::validators::get_single_string(values);

        if (s == phase_t::PHASE_BUILD || s == phase_t::PHASE_PLACE)
        {
            v = boost::any(phase_t(s));
        }
        else
        {
            throw po::validation_error(po::validation_error::invalid_option_value);
        }
    }

    struct model_t
    {
        static std::set<std::string> values;
        model_t(const std::string& v) : value(v) { }
        std::string value;
    };


    void validate(boost::any& v, std::vector<std::string> const& values, model_t*, int)
    {
        // Make sure no previous assignment to 'v' was made.
        po::validators::check_first_occurrence(v);

        // Extract the first string from 'values'. If there is more than
        // one string, it's an error, and exception will be thrown.
        std::string const& s = po::validators::get_single_string(values);

        if (ar::evo_models.find(s) != ar::evo_models.end())
        {
            v = boost::any(model_t(s));
        }
        else
        {
            throw po::validation_error(po::validation_error::invalid_option_value);
        }
    }

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
                ((PHASE + "," + PHASE_SHORT).c_str(), po::value<phase_t>()->required(),
                        "One of 'b' for 'Build' or 'p' for 'Place'"
                        "* b: Build DB of phylo-kmers (done 1 time)."
                        "* p: Phylogenetic placement itself (done n times)"
                        "requires the DB generated during 'build' phase.")
                ((ARBINARY + "," + ARBINARY_SHORT).c_str(), po::value<fs::path>()->required(),
                       "[file] Binary for marginal AR, currently 'phyml' and 'baseml' (from PAML) are supported. (b phase)")
                ((REFALIGN + "," + REFALIGN_SHORT).c_str(), po::value<fs::path>()->required(),
                        "Reference alignment file")
                ((REFTREE + "," + REFTREE_SHORT).c_str(), po::value<fs::path>()->required(),
                        "Phylogenetic tree file")
                ((WORKING_DIR + "," + WORKING_DIR_SHORT).c_str(), po::value<fs::path>()->default_value(fs::current_path()),
                        "Path to the working directory (b|p phase).")
                ((ARDIR + "," + ARDIR_SHORT).c_str(), po::value<fs::path>(),
                        "Ancestral reconstruction results")
                ((K + "," + K_SHORT).c_str(), po::value<size_t>()->default_value(8),
                         "k-mer length used at DB build. (b phase)")
                ((MODEL + "," + MODEL_SHORT).c_str(), po::value<model_t>()->required(),
                        "Model used in AR, one of the following:\n\t*nucl: JC69, HKY85, K80, F81, TN93, GTR")
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
            else if (vm.count(PHASE))
            {
                phase_t phase = vm[PHASE].as<phase_t>();
                if (phase.value == phase_t::PHASE_BUILD)
                {
                    parameters.action = build;
                }
                else if (phase.value == phase_t::PHASE_PLACE)
                {
                    parameters.action = place;
                }
            }

            if (vm.count(ARBINARY))
            {
                parameters.ar_binary = vm[ARBINARY].as<fs::path>().string();
            }
            if (vm.count(REFTREE))
            {
                parameters.tree_file = vm[REFTREE].as<fs::path>().string();
            }
            if (vm.count(REFALIGN))
            {
                parameters.alignment_file = vm[REFALIGN].as<fs::path>().string();
            }
            parameters.working_directory = vm[WORKING_DIR].as<fs::path>().string();

            if (vm.count(ARDIR))
            {
                parameters.ar_directory = vm[ARDIR].as<fs::path>().string();
            }

            model_t model = vm[MODEL].as<model_t>();
            const float alpha = 1.0f;
            const int categories = 4;
            parameters.model = ar::evo_model(model.value, alpha, categories);
            parameters.kmer_size = vm[K].as<size_t>();
            parameters.ar_parameters = "";
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

