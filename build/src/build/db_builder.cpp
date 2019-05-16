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
using namespace core;

namespace rappas
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(const string& working_directory, const string& ar_probabilities_file,
                                   const string& original_tree_file, const string& extended_tree_file,
                                   const string& extended_mapping_file, const string& artree_mapping_file, size_t kmer_size);
    public:
        /// Member types
        /// \brief A hash map to store all the phylo-kmers, placed to one original node
        using branch_hash_map = hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;

        /// \brief A group of node ids that must be processed together. We group together
        /// the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(const string& working_directory, const string& ar_probabilities_file,
                   const string& original_tree_file, const string& extended_tree_file,
                   const string& extended_mapping_file, const string& artree_mapping_file, size_t kmer_size);
        db_builder(const db_builder&) = delete;
        db_builder(db_builder&&) = delete;
        db_builder& operator=(const db_builder&) = delete;
        db_builder& operator=(db_builder&&) = delete;
        ~db_builder() noexcept = default;


        /// \brief Runs the database construction
        void run();

    private:
        /// \brief Groups ghost nodes by corresponding original node id
        std::vector<id_group> group_ghost_ids(const std::vector<std::string>& ghost_ids) const;

        /// \brief Groups references to submatrices of probabilites, corresponding to a group of nodes
        proba_group get_submatrices(const proba_matrix& probas, const id_group& group) const;

        /// \brief Runs a phylo-kmer exploration for every ghost node of the extended_tree
        /// \return The number of explored phylo-kmers. This number can be more, than a size of a resulting database
        size_t explore_kmers(const phylo_tree& orignal_tree, const phylo_tree& extended_tree, const proba_matrix& probas);

        /// \brief Explores phylo-kmers of a collection of ghost nodes. Here we assume that the nodes
        /// in the group correspond to one original node
        /// \return A hash map with phylo-kmers stored and a number of explored phylo-kmers
        std::pair<branch_hash_map, size_t> explore_group(const proba_group& group);

        string _working_directory;
        string _ar_probabilities_file;
        string _original_tree_file;
        string _extended_tree_file;
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

db_builder::db_builder(const string& working_directory, const string& ar_probabilities_file, const string& original_tree_file,
    const string& extended_tree_file, const string& extended_mapping_file, const string& artree_mapping_file, size_t kmer_size)
    : _working_directory{ working_directory }
    , _ar_probabilities_file{ ar_probabilities_file }
    , _original_tree_file{ original_tree_file }
    , _extended_tree_file{ extended_tree_file }
    , _extended_mapping_file{ extended_mapping_file }
    , _artree_mapping_file{ artree_mapping_file }
    , _kmer_size{ kmer_size }
    , _phylo_kmer_db{ kmer_size }
{}

void db_builder::run()
{
    /// Load .tsv files
    _extended_mapping = rappas::io::load_extended_mapping(_extended_mapping_file);
    _artree_mapping = rappas::io::load_artree_mapping(_artree_mapping_file);

    /// Load .newick files
    const auto original_tree = rappas::io::load_newick(_original_tree_file);
    const auto extended_tree = rappas::io::load_newick(_extended_tree_file);

    /// Load PhyML output
    const auto proba_matrix = rappas::io::load_phyml_probas(_ar_probabilities_file);

    /// Run the branch and bound algorithm
    std::cout << "Building database..." << std::endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto tuples_count = explore_kmers(original_tree, extended_tree, proba_matrix);
    auto end = std::chrono::steady_clock::now();

    /// Calculate the number of phylo-kmers stored in the database
    size_t total_entries = 0;
    for (const auto& kmer_entry : _phylo_kmer_db)
    {
        total_entries += kmer_entry.second.size();
    }

    std::cout << "Built " << total_entries << " phylo-kmers out of " << tuples_count << " for " << _phylo_kmer_db.size()
              << " k-mer values.\nTime (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "\n\n" << std::flush;
}

bool is_ghost(const phylo_node& node)
{
    const string& label = node.get_label();
    return boost::ends_with(label, "_X0") || boost::ends_with(label, "_X1");
}

/// \brief Returns a list of ghost node ids
std::vector<std::string> get_ghost_ids(const phylo_tree& tree)
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

std::vector<db_builder::id_group> db_builder::group_ghost_ids(const std::vector<std::string>& ghost_ids) const
{
    std::vector<id_group> groups;
    groups.reserve(ghost_ids.size() / 2);

    std::unordered_map<branch_type, size_t> mapping;
    for (const auto& ghost_id : ghost_ids)
    {
        const auto& original_preorder_id = _extended_mapping.at(ghost_id);
        if (const auto it = mapping.find(original_preorder_id); it != mapping.end())
        {
            groups[it->second].push_back(ghost_id);
        }
        else
        {
            groups.push_back({ ghost_id });
            mapping[original_preorder_id] = groups.size() - 1;
        }

    }
    return groups;
}

db_builder::proba_group db_builder::get_submatrices(const proba_matrix& probas, const id_group& group) const
{
    proba_group submatrices;

    for (const auto& branch_node_label : group)
    {
        const auto& artree_node_label = _artree_mapping.at(branch_node_label);
        if (const auto& it = probas.find(artree_node_label); it != probas.end())
        {
            submatrices.push_back(std::cref(it->second));
        }
        else
        {
            std::cerr << "Internal error: could not find " << artree_node_label << " node." << std::endl;
        }
    }

    return submatrices;
}

size_t db_builder::explore_kmers(const phylo_tree& orignal_tree, const phylo_tree& extended_tree, const proba_matrix& probas)
{
    size_t count = 0;

    /// Here we assume that every original node corresponds to two ghost nodes.
    const size_t ghosts_per_node = 2;

    /// Filter and group ghost nodes
    const auto node_groups = group_ghost_ids(get_ghost_ids(extended_tree));

    /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
    /// in a hash map for every group separately.
    _branch_maps.resize(orignal_tree.get_node_count());
    std::vector<phylo_kmer::branch_type> node_postorder_ids(orignal_tree.get_node_count());

    //#pragma omp parallel for schedule(auto) reduction(+: count) //num_threads(8)
    for (size_t i = 0; i < node_groups.size(); ++i)
    {
        const auto& node_group = node_groups[i];
        assert(node_group.size() == ghosts_per_node);
        (void)ghosts_per_node;

        /// Having a label of node in the extended tree, we need to find the corresponding node
        /// in the original tree. We take the first ghost node, because all of them correspond the same
        /// original node
        const auto original_node_preorder_id = _extended_mapping[node_group[0]];
        const auto original_node_postorder_id = (*orignal_tree[original_node_preorder_id])->get_postorder_id();
        node_postorder_ids[i] = original_node_postorder_id;

        /// Get submatrices of probabilities for a group
        proba_group submatrices = get_submatrices(probas, node_group);
        size_t branch_count = 0;

        /// Explore k-mers of a group and save results in a hash map
        std::tie(_branch_maps[i], branch_count) = explore_group(submatrices);
        count += branch_count;
    }

    /// Merge hash maps in a final data structure
    for (size_t i = 0; i < _branch_maps.size(); ++i)
    {
        auto& map = _branch_maps[i];
        for (const auto& [key, score] : map)
        {
            _phylo_kmer_db.insert(key, { node_postorder_ids[i], score });
        }
        /// Replace a map with an empty one to free memory
        map = {};
    }

    return count;
}

/// \brief Puts a key-value pair in a hash map. Used to process branches in parallel
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

std::pair<db_builder::branch_hash_map, size_t> db_builder::explore_group(const proba_group& group)
{
    branch_hash_map hash_map;
    size_t count = 0;

    for (auto node_entry_ref : group)
    {
        const auto& node_entry = node_entry_ref.get();
        for (auto window = node_entry.begin(_kmer_size); window != node_entry.end(); ++window)
        {
            for (const auto& kmer : *window)
            {
                put(hash_map, kmer.key, kmer.score);
                ++count;
            }
        }
    }

    return { std::move(hash_map), count };
}

namespace rappas
{
    phylo_kmer_db build(const std::string& working_directory, const std::string& ar_probabilities_file,
                        const std::string& original_tree_file, const std::string& extended_tree_file,
                        const std::string& extended_mapping_file, const std::string& artree_mapping_file, size_t kmer_size)
    {
        db_builder builder(working_directory, ar_probabilities_file, original_tree_file, extended_tree_file,
                           extended_mapping_file, artree_mapping_file, kmer_size);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
