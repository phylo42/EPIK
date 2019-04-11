#include "db_builder.h"
#include "tree/phylo_tree.h"
#include "pp_matrix/proba_matrix.h"
#include "pp_matrix/phyml.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/log/trivial.hpp>
#include <chrono>
#include <iomanip>

namespace fs = boost::filesystem;
using std::string;
using std::vector;
using std::cout, std::endl;
using std::to_string;

db_builder::db_builder(string working_directory, string ar_probabilities_file, string tree_file,
    string extended_mapping_file, string artree_mapping_file, size_t kmer_size)
    : _working_directory(move(working_directory))
    , _ar_probabilities_file(move(ar_probabilities_file))
    , _tree_file(move(tree_file))
    , _extended_mapping_file(move(extended_mapping_file))
    , _artree_mapping_file(move(artree_mapping_file))
    , _kmer_size(kmer_size)
{}

size_t db_builder::explore_branch(const branch_entry& probas, branch_node_t original_id)
{
    size_t count = 0;
    for (auto window = probas.begin(_kmer_size); window != probas.end(); ++window)
    {
        for (const auto& kmer : *window)
        {
            _phylo_kmer_db.put(kmer.value, original_id, kmer.score);
            ++count;
        }
    }
    return count;
}

void db_builder::explore_kmers(const phylo_tree& tree, const proba_matrix& probas)
{
    size_t count = 0;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    /// iterate over fake nodes
    for (const auto& branch_node: tree)
    {
        if (is_fake(branch_node))
        {
            const auto original_id = _extended_mapping[branch_node.get_label()];

            /// get submatrix of probabilities for a current branch node (if presented in proba matrix)
            const auto phyml_branch_label = _artree_mapping[branch_node.get_label()];
            if (const auto& it = probas.find(phyml_branch_label); it != probas.end())
            {
                count += explore_branch(it->second, original_id);
            }
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "Phylokmer generation time (s) = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
    std::cout << "Phylokmer generation time (ms) = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;

    size_t total_entries = 0;
    for (const auto& kmer_entry : _phylo_kmer_db)
    {
        total_entries += kmer_entry.second.size();
    }
    std::cout << "Kmers generated: " << _phylo_kmer_db.size() << "\n";
    std::cout << "Kmer tuples in the hash: " << total_entries << "\n";
    std::cout << "Tuples explored: " << count << "\n";
}

return_code db_builder::run()
{
    _extended_mapping = load_extended_mapping(_extended_mapping_file);
    _artree_mapping = load_artree_mapping(_artree_mapping_file);

    const auto tree = load_newick(_tree_file);
    const auto probas = load_phyml_probas(_ar_probabilities_file);
    explore_kmers(tree, probas);

    return return_code::success;
}