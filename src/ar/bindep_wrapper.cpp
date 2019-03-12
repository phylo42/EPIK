#include "bindep_wrapper.h"
#include "evo_model.h"
#include <vector>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
#include <iostream>

namespace fs = boost::filesystem;

namespace ar
{
    bindep_wrapper::bindep_wrapper(const std::string &ar_binary, const std::string &alignment_file,
                                   const std::string &tree_file, const std::string &ar_parameters,
                                   const ar::evo_model &model)
        : command(build_phyml_command(ar_binary, alignment_file, tree_file, ar_parameters, model))
    { }

    bindep_wrapper::~bindep_wrapper() { }

    void bindep_wrapper::run() const
    {
        std::cout << command << std::endl;
    }

    std::string bindep_wrapper::build_phyml_command(const std::string &ar_binary, const std::string &alignment_file,
                                                    const std::string &tree_file, const std::string &ar_parameters,
                                                    const ar::evo_model &model) const
    {
        std::vector<std::string> command =
                {
                        // binary file itself
                        fs::canonical(ar_binary, fs::current_path()).string(),

                        // marginal reconstruct
                        "--ancestral",

                        // no interactive questions for mem usage
                        "--no_memory_check",

                        // alignment file
                        "-i", fs::canonical(alignment_file).string(),

                        // tree file
                        "-u", fs::canonical(tree_file).string()
                };

        if (ar_parameters.empty())
        {
            // model
            command.emplace_back("-m");
            command.push_back(model.name);
            /*
            if (model.isProteinModel()) {
                com.add("-d"); //analysis type
                com.add("aa");
            }*/

            // number of relative substitution rate categories
            command.emplace_back("-c");
            command.push_back(std::to_string(model.categories));

            // neither approximate likelihood ratio test nor bootstrap values are computed
            command.emplace_back("-b");
            command.emplace_back("0");

            // proportion of invariable sites
            command.emplace_back("-v");
            command.emplace_back("0.0");

            // rate parameters are optimised
            command.emplace_back("-o");
            command.emplace_back("r");

            // gamma shape param
            command.emplace_back("-a");
            command.push_back(std::to_string(model.alpha));

            // base frequencies based on aligned
            command.emplace_back("-f");
            command.emplace_back("e");

            // no interactive questions
            // command.emplace_back("--quiet");

            // phyml version from 06.2018 do not accept anymore duplicate seqs
            // TODO: find a way to test and choose correct parameters depending on version
            // command.emplace_back("--leave_duplicates");
        }
        else
        {
            throw std::runtime_error("ar_parameters are not supported yet");
            // if parameters given by user via --arparameters, forget previous
            // command and use this one instead.
            // com.addAll(Arrays.asList(ARParameters.split(" ")));
        }
        return boost::algorithm::join(command, " ");
    }
}