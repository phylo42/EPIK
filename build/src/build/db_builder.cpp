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
using core::phylo_kmer, core::phylo_kmer_db, core::phylo_tree;

namespace rappas
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(const string& working_directory, const string& ar_probabilities_file,
                                   const string& tree_file, const string& extended_mapping_file,
                                   const string& artree_mapping_file, size_t kmer_size);
    public:
        using branch_hash_map = core::hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;

        db_builder(const string& working_directory, const string& ar_probabilities_file,
                   const string& tree_file, const string& extended_mapping_file,
                   const string& artree_mapping_file, size_t kmer_size);

        void run();

    private:
        size_t explore_kmers(const phylo_tree& tree, const proba_matrix& probas);
        std::pair<branch_hash_map, size_t> explore_branch(const node_entry& probas);

        string _working_directory;
        string _ar_probabilities_file;
        string _tree_file;
        string _extended_mapping_file;
        string _artree_mapping_file;

        size_t _kmer_size;
        phylo_kmer_db _phylo_kmer_db;
        std::vector<branch_hash_map> _branch_maps;

        extended_mapping _extended_mapping;
        artree_label_mapping _artree_mapping;
    };

}

using namespace rappas;

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

/// Puts a key-value pair in a hash map. Used to process branches in parallel
void put(db_builder::branch_hash_map& map, phylo_kmer::key_type key, phylo_kmer::score_type score)
{
    if (auto it = map.find(key); it != map.end())
    {
        if (it->second < score)
        {
            map[key] = score;
        }
    }
    else
    {
        map[key] = score;
    }
}

std::pair<db_builder::branch_hash_map, size_t> db_builder::explore_branch(const node_entry& probas)
{
    branch_hash_map hash_map;
    size_t count = 0;
    for (auto window = probas.begin(_kmer_size); window != probas.end(); ++window)
    {
        for (const auto& kmer : *window)
        {
            put(hash_map, kmer.key, kmer.score);
            ++count;
        }
    }
    return { std::move(hash_map), count };
}

bool is_ghost(const core::phylo_node& node)
{
    const string& label = node.get_label();
    return boost::ends_with(label, "_X0") || boost::ends_with(label, "_X1");
}

std::vector<std::string> get_ghost_ids(const core::phylo_tree& tree)
{
    std::vector<std::string> branch_ids;

    for (const auto& branch_node: tree)
    {
        if (is_ghost(branch_node))
        {
            branch_ids.push_back(branch_node.get_label());
        }
    }
    return branch_ids;
}


size_t db_builder::explore_kmers(const core::phylo_tree& tree, const proba_matrix& probas)
{
    size_t count = 0;

    /// Filter ghost nodes
    const auto ghost_node_ids = get_ghost_ids(tree);
    std::vector<phylo_kmer::branch_type> original_node_ids(ghost_node_ids.size());

    /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
    /// in a hash map for every branch separately.
    _branch_maps.resize(tree.get_node_count());
    #pragma omp parallel for schedule(auto) // num_threads(num_threads)
    for (size_t i = 0; i < ghost_node_ids.size(); ++i)
    {
        const auto branch_node_label = ghost_node_ids[i];
        original_node_ids[i] = _extended_mapping[branch_node_label];

        /// Get submatrix of probabilities for a current branch node (if presented in proba matrix)
        const auto phyml_node_label = _artree_mapping[branch_node_label];
        if (const auto& it = probas.find(phyml_node_label); it != probas.end())
        {
            size_t branch_count;
            std::tie(_branch_maps[i], branch_count) = explore_branch(it->second);
            count += branch_count;
        }
    }

    /// Merge hash maps in a final data structure
    for (size_t i = 0; i < _branch_maps.size(); ++i)
    {
        auto& map = _branch_maps[i];
        for (const auto& [key, score] : map)
        {
            _phylo_kmer_db.put(key, original_node_ids[i], score);
        }
        /// Replace a map with an empty one to free memory
        map = {};
    }

    return count;
}

void db_builder::run()
{
    _extended_mapping = rappas::io::load_extended_mapping(_extended_mapping_file);
    _artree_mapping = rappas::io::load_artree_mapping(_artree_mapping_file);

    const auto tree = rappas::io::load_newick(_tree_file);
    const auto proba_matrix = rappas::io::load_phyml_probas(_ar_probabilities_file);

    /// Run the branch and bound algorithm
    std::cout << "Building database..." << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto tuples_count = explore_kmers(tree, proba_matrix);
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

namespace rappas
{
    phylo_kmer_db build(const std::string& working_directory, const std::string& ar_probabilities_file,
                        const std::string& tree_file, const std::string& extended_mapping_file,
                        const std::string& artree_mapping_file, size_t kmer_size)
    {
        db_builder builder(working_directory, ar_probabilities_file, tree_file,
                           extended_mapping_file, artree_mapping_file, kmer_size);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
