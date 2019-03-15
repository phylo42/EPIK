#include "db_builder.h"
#include "phylo_tree.h"
#include "ar.h"
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
           size_t kmer_size,
           const seq_traits& traits)
    : _working_directory(working_directory)
    , _ar_probabilities_file(ar_probabilities_file)
    , _tree_file(tree_file)
    , _kmer_size(kmer_size)
    , _seq_traits(traits)
{}

return_code_t db_builder::run()
{
    phylo_tree tree = load_newick(_tree_file);
    proba_matrix probas = load_phyml_probas(_ar_probabilities_file);
    return return_code::success;
}