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

namespace fs = boost::filesystem;
using std::string;
using std::vector;
using std::cout, std::endl;
using std::to_string;


/// Apply a node mapping to a phylogenetic tree
void apply_mapping(phylo_tree& tree, const node_mapping& mapping)
{
    for (auto& node : tree)
    {
        auto it = mapping.from_phyml.find(node.get_label());
        if (it != end(mapping.from_phyml))
        {
            node.set_label(it->second);
        }
    }
}

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

void db_builder::explore_branch(const branch_entry& branch)
{
    std::cout << "Exploring branch " << branch.get_branch_id() << std::endl;

    for (auto window = branch.begin(_kmer_size); window != branch.end(); ++window)
    {
        for (auto kmer : *window)
        {
            std::cout << kmer.value << "\t" << kmer.score << "\t" << branch.get_branch_id()
                << '\t' << window->get_start_pos() << '\n';
        }
    }
}

return_code_t db_builder::run()
{
    const node_mapping mapping = load_node_mapping(_mapping_file);
    auto tree = load_newick(_tree_file);

    /// restore original names of inner branch nodes that has been rewritten by PhyML
    /// TODO: encapsulate this in a AR strategy class
    apply_mapping(tree, mapping);

    const auto probas = load_phyml_probas(_ar_probabilities_file, _seq_traits);

    /// iterate over fake nodes
    for (const auto& branch_node: std::as_const(tree))
    {
        if (is_fake(branch_node))
        {
            /// get submatrix of probabilities for a current branch node
            auto phyml_branch_id = std::stoi(mapping.to_phyml.at(branch_node.get_label()));
            explore_branch(probas.at(phyml_branch_id));
        }
    }
    return return_code::success;
}