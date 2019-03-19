#include "db_builder.h"
#include "phylo_tree.h"
#include "ar.h"
#include "phylo_kmer_explorer.h"
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

return_code_t db_builder::run()
{
    phylo_tree tree = load_newick(_tree_file);
    proba_matrix probas = load_phyml_probas(_ar_probabilities_file);
    node_mapping mapping = load_node_mapping(_mapping_file);

    for (const auto& branch_nodes: tree)
    {
        if (is_fake(branch_nodes))
        {
            cout << branch_nodes.get_label() << endl;
        }

        /*
        probas_view view = probas.make_view(branch_node);

        phylo_kmer_explorer kmer_explorer(_seq_traits, view, branch_node);
        while (kmer_explorer.has_next())
        {
            phylo_kmer ph_kmer = kmer_explorer.next_phylo_kmer();
            _phylo_kmer_db[ph_kmer.kmer_value] = ph_kmer;
        }
*/
    }
    return_code::success;
}