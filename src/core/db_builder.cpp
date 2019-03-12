#include "db_builder.h"
#include "alignment.h"
#include "newick.h"
#include "fasta.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>

namespace fs = boost::filesystem;
using std::string;
using std::vector;
using std::cout, std::endl;
using std::to_string;

db_builder::db_builder(const string& ar_binary, const string& alignment_file,
                       const string& tree_file, const string& working_directory,
                       const string& ar_directory, const string& ar_parameters,
                       size_t kmer_size, const ar::evo_model& model,
                       const seq_traits& traits, double reduction_ratio)
    : _ar_binary(ar_binary)
    , _alignment_file(alignment_file)
    , _tree_file(tree_file)
    , _working_directory(working_directory)
    , _ar_directory(ar_directory)
    , _ar_parameters(ar_parameters)
    , _kmer_size(kmer_size)
    , _model(model)
    , _seq_traits(traits)
    , _reduction_ratio(reduction_ratio)
{}

return_code_t db_builder::run()
{
    //std::system(build_phyml_command().c_str());

    // read and validate the alignment
    alignment align = load_alignment(_alignment_file);
    alignment_validation validator(_seq_traits);
    validator.validate(align);

    // reduce_alignment the alignment and save it to the working directory
    alignment reduced_alignment = reduce_alignment(_seq_traits, align, _reduction_ratio);
    string alignment_output_file = fs::weakly_canonical(fs::path(_working_directory) / fs::path("align.reduced")).string();
    save_alignment(reduced_alignment, alignment_output_file);

    cout << "Alignment reduction based on ratio = " << _reduction_ratio << endl;
    cout << "Before reduction: Dimension = " << align.width() << "x" << align.height() << endl;
    cout << "After reduction: Dimension = " << reduced_alignment.width() << "x" << reduced_alignment.height() << endl;


    phylo_tree tree = load_newick(_tree_file);
    cout << tree.get_node_count() << endl;

    //cout << build_phyml_command() << endl;
    return return_code::success;
}

string db_builder::build_phyml_command() const
{
    vector<string> command =
    {
        fs::canonical(_ar_binary, fs::current_path()).string(),
        "--ancestral", // marginal reconstruct
        "--no_memory_check", // no interactive questions for mem usage
        "-i", // alignment file
        fs::canonical(_alignment_file).string(),
        "-u", // tree file
        fs::canonical(_tree_file).string()
    };

    if (_ar_parameters.empty()) {
        // model
        command.emplace_back("-m");
        command.push_back(_model.name);
        /*
        if (model.isProteinModel()) {
            com.add("-d"); //analysis type
            com.add("aa");
        }*/

        // number of relative substitution rate categories
        command.emplace_back("-c");
        command.push_back(to_string(_model.categories));

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
        command.push_back(to_string(_model.alpha));

        // base frequencies based on aligned
        command.emplace_back("-f");
        command.emplace_back("e");

        // no interactive questions
        // command.emplace_back("--quiet");

        // phyml version from 06.2018 do not accept anymore duplicate seqs
        // TODO: find a way to test and choose correct parameters depending on version
        // command.emplace_back("--leave_duplicates");
    } else {
        throw std::runtime_error("ar_parameters are not supported yet");
        // if parameters given by user via --arparameters, forget previous
        // command and use this one instead.
        // com.addAll(Arrays.asList(ARParameters.split(" ")));
    }
    return boost::algorithm::join(command, " ");
}
