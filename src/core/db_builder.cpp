#include "db_builder.h"
#include "core/tree/phylo_tree.h"
#include "pp_matrix/proba_matrix.h"
#include "core/pp_matrix/phyml.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <chrono>

namespace fs = boost::filesystem;
using std::string;
using std::vector;
using std::cout, std::endl;
using std::to_string;

db_builder::db_builder(const std::string& working_directory,
           const std::string& ar_probabilities_file,
           const std::string& tree_file,
           const std::string& mapping_file,
           size_t kmer_size,
           const seq_traits& traits)
    : _working_directory(working_directory)
    , _ar_probabilities_file(ar_probabilities_file)
    , _tree_file(tree_file)
    , _mapping_file(mapping_file)
    , _kmer_size(kmer_size)
    , _seq_traits(traits)
{}

void db_builder::explore_branch(const branch_entry& probas)
{
    const auto artree_branch_id = _mapping.artree_label_to_artree_id[std::to_string(probas.get_branch_label())];
    //std::cout << "Exploring branch " << artree_branch_id << " (artree_label " << probas.get_branch_label() << ")\n";
    for (auto window = probas.begin(_kmer_size); window != probas.end(); ++window)
    {
        for (auto kmer : *window)
        {
            /*std::cout << kmer.value << "\t" << kmer.score << "\t" << artree_branch_id
                << '\t' << window->get_start_pos() << '\n' << std::flush;*/
            (void)kmer;
        }
    }
}

return_code_t db_builder::run()
{
    _mapping = load_node_mapping(_mapping_file);
    auto tree = load_newick(_tree_file);

    //apply_mapping(tree, mapping);

    const auto probas = load_phyml_probas(_ar_probabilities_file, _seq_traits);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    /// iterate over fake nodes
    for (const auto& branch_node: std::as_const(tree))
    {
        if (is_fake(branch_node))
        {
            /// get submatrix of probabilities for a current branch node (if presented in proba matrix)
            const auto phyml_branch_label = branch_id( std::stoul(_mapping.extended_label_to_phyml_label.at(branch_node.get_label())));
            if (auto it = probas.find(phyml_branch_label); it != probas.end())
            {
                explore_branch(it->second);
            }
        }
    }
    std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    std::cout << "Phylokmer generation time (s) = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<std::endl;
    std::cout << "Phylokmer generation time (ms) = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;
    return return_code::success;
}