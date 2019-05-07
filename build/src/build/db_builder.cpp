#include <iostream>
#include <chrono>
#include <boost/algorithm/string/predicate.hpp>

#include <core/phylo_kmer_db.h>
#include <core/phylo_tree.h>
#include "db_builder.h"
#include "pp_matrix/proba_matrix.h"
#include "pp_matrix/phyml.h"

using std::string;
using std::cout, std::endl;
using std::to_string;

class db_builder
{
    friend core::phylo_kmer_db build_database(const std::string& working_directory, const std::string& ar_probabilities_file,
                                              const std::string& tree_file, const std::string& extended_mapping_file,
                                              const std::string& artree_mapping_file, size_t kmer_size);
public:
    db_builder(const std::string& working_directory, const std::string& ar_probabilities_file,
               const std::string& tree_file, const std::string& extended_mapping_file,
               const std::string& artree_mapping_file, size_t kmer_size);

    void run();

private:
    size_t explore_kmers(const core::phylo_tree& tree, const proba_matrix& probas);
    size_t explore_branch(const branch_entry& probas, core::phylo_kmer::branch_type common_branch_label);

    std::string _working_directory;
    std::string _ar_probabilities_file;
    std::string _tree_file;
    std::string _extended_mapping_file;
    std::string _artree_mapping_file;

    size_t _kmer_size;
    core::phylo_kmer_db _phylo_kmer_db;
    extended_mapping _extended_mapping;
    artree_label_mapping _artree_mapping;
};



db_builder::db_builder(const string& working_directory, const string& ar_probabilities_file, const string& tree_file,
                       const string& extended_mapping_file, const string& artree_mapping_file, size_t kmer_size)
    : _working_directory{ working_directory }
    , _ar_probabilities_file{ ar_probabilities_file }
    , _tree_file{ tree_file }
    , _extended_mapping_file{ extended_mapping_file }
    , _artree_mapping_file{ artree_mapping_file }
    , _kmer_size{ kmer_size }
    , _phylo_kmer_db{ kmer_size }
{}

size_t db_builder::explore_branch(const branch_entry& probas, core::phylo_kmer::branch_type original_id)
{
    size_t count = 0;
    for (auto window = probas.begin(_kmer_size); window != probas.end(); ++window)
    {
        for (const auto& kmer : *window)
        {
            _phylo_kmer_db.put(kmer.key, original_id, kmer.score);
            ++count;
        }
    }
    return count;
}

bool is_fake(const core::phylo_node& node)
{
    const string& label = node.get_label();
    return boost::ends_with(label, "_X0") || boost::ends_with(label, "_X1");
}

size_t db_builder::explore_kmers(const core::phylo_tree& tree, const proba_matrix& probas)
{
    size_t count = 0;

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
    return count;
}

void db_builder::run()
{
    _extended_mapping = load_extended_mapping(_extended_mapping_file);
    _artree_mapping = load_artree_mapping(_artree_mapping_file);

    const auto tree = core::load_newick(_tree_file);
    const auto probas = load_phyml_probas(_ar_probabilities_file);

    /// Run the branch and bound algorithm
    std::cout << "Building database..." << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto tuples_count = explore_kmers(tree, probas);
    auto end = std::chrono::steady_clock::now();

    size_t total_entries = 0;
    for (const auto& kmer_entry : _phylo_kmer_db)
    {
        total_entries += kmer_entry.second.size();
    }

    std::cout << "Built " << total_entries << " phylo-kmers out of " << tuples_count << " for " << _phylo_kmer_db.size()
        << " k-mer values.\nTime (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
        << "\n\n" << std::flush;
}

core::phylo_kmer_db build_database(const std::string& working_directory, const std::string& ar_probabilities_file,
                                   const std::string& tree_file, const std::string& extended_mapping_file,
                                   const std::string& artree_mapping_file, size_t kmer_size)
{
    db_builder builder(working_directory, ar_probabilities_file, tree_file,
        extended_mapping_file, artree_mapping_file, kmer_size);
    builder.run();
    return std::move(builder._phylo_kmer_db);
}