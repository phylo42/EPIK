#include <iostream>
#include <chrono>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <core/phylo_kmer_db.h>
#include <core/phylo_tree.h>
#include <utils/io/file_io.h>
#include "db_builder.h"
#include "pp_matrix/proba_matrix.h"
#include "pp_matrix/phyml.h"

using std::string;
using std::cout, std::endl;
using std::to_string;
using namespace core;
namespace fs = boost::filesystem;


namespace rappas
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(const string& working_directory, const string& ar_probabilities_file,
                                   const string& original_tree_file, const string& extended_tree_file,
                                   const string& extended_mapping_file, const string& artree_mapping_file,
                                   size_t kmer_size, core::phylo_kmer::score_type omega, size_t num_threads);
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
                   const string& extended_mapping_file, const string& artree_mapping_file,
                   size_t kmer_size, core::phylo_kmer::score_type omega, size_t num_threads);
        db_builder(const db_builder&) = delete;
        db_builder(db_builder&&) = delete;
        db_builder& operator=(const db_builder&) = delete;
        db_builder& operator=(db_builder&&) = delete;
        ~db_builder() noexcept = default;


        /// \brief Runs the database construction
        void run();

    private:
        /// \brief Returns a filename for a hashmap of a given group
        std::string group_hashmap_file(const branch_type& group) const;

        /// \brief Groups ghost nodes by corresponding original node id
        std::vector<id_group> group_ghost_ids(const std::vector<std::string>& ghost_ids) const;

        /// \brief Groups references to submatrices of probabilites, corresponding to a group of nodes
        proba_group get_submatrices(const proba_matrix& probas, const id_group& group) const;

        /// \brief Runs a phylo-kmer exploration for every ghost node of the extended_tree
        /// \return The number of explored phylo-kmers. This number can be more, than a size of a resulting database
        size_t explore_kmers(const phylo_tree& original_tree, const phylo_tree& extended_tree,
                             const proba_matrix& probas);

        /// \brief Explores phylo-kmers of a collection of ghost nodes. Here we assume that the nodes
        /// in the group correspond to one original node
        /// \return A hash map with phylo-kmers stored and a number of explored phylo-kmers
        std::pair<branch_hash_map, size_t> explore_group(const proba_group& group);

        /// \brief Saves a hash map to file
        void save_hash_map(const branch_hash_map& map, const std::string& filename) const;

        /// \brief Loads a hash map from file
        branch_hash_map load_hash_map(const std::string& filename) const;

        /// \brief Working and output directory
        string _working_directory;
        /// \brief A subdirectory of _working_directory to store hashmaps
        string _hashmaps_directory;

        string _ar_probabilities_file;
        string _original_tree_file;
        string _extended_tree_file;
        string _extended_mapping_file;
        string _artree_mapping_file;

        size_t _kmer_size;
        core::phylo_kmer::score_type _omega;
        size_t _num_threads;
        phylo_kmer_db _phylo_kmer_db;

        extended_mapping _extended_mapping;
        artree_label_mapping _artree_mapping;
    };
}

using namespace rappas;

db_builder::db_builder(const string& working_directory, const string& ar_probabilities_file,
    const string& original_tree_file, const string& extended_tree_file, const string& extended_mapping_file,
    const string& artree_mapping_file, size_t kmer_size, core::phylo_kmer::score_type omega, size_t num_threads)
    : _working_directory{ working_directory }
    , _hashmaps_directory{ (fs::path{ working_directory } / fs::path{ "hashmaps" }).string() }
    , _ar_probabilities_file{ ar_probabilities_file }
    , _original_tree_file{ original_tree_file }
    , _extended_tree_file{ extended_tree_file }
    , _extended_mapping_file{ extended_mapping_file }
    , _artree_mapping_file{ artree_mapping_file }
    , _kmer_size{ kmer_size }
    , _omega{ omega }
    , _num_threads{ num_threads }
    /// I do not like reading a file here, but it seems to be better than having something like "set_tree"
    /// in the public interface of the phylo_kmer_db class.
    , _phylo_kmer_db{ kmer_size, omega, rappas::io::read_as_string(original_tree_file) }
{}

void create_directory(const std::string& dirname)
{
    /// Create if does not exist
    if (!fs::is_directory(dirname) || !fs::exists(dirname))
    {
        if (!fs::create_directory(dirname))
        {
            throw new std::runtime_error("Cannot create directory " + dirname);
        }
    }
}

void db_builder::run()
{
    /// Load .tsv files
    _extended_mapping = rappas::io::load_extended_mapping(_extended_mapping_file);
    _artree_mapping = rappas::io::load_artree_mapping(_artree_mapping_file);

    /// create a temporary directory for hashmaps
    create_directory(_hashmaps_directory);

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

    std::cout << "Built " << total_entries << " phylo-kmers out of " << tuples_count << " for "
              << _phylo_kmer_db.size() << " k-mer values.\nTime (ms): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
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

std::string db_builder::group_hashmap_file(const branch_type& group) const
{
    return (fs::path{ _hashmaps_directory } / fs::path{ std::to_string(group) }).string();
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

size_t db_builder::explore_kmers(const phylo_tree& original_tree, const phylo_tree& extended_tree,
                                 const proba_matrix& probas)
{
    size_t count = 0;

    /// Here we assume that every original node corresponds to two ghost nodes.
    const size_t ghosts_per_node = 2;

    /// Filter and group ghost nodes
    const auto node_groups = group_ghost_ids(get_ghost_ids(extended_tree));

    /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
    /// in a hash map for every group separately on disk.
    std::vector<phylo_kmer::branch_type> node_postorder_ids(node_groups.size());

    #pragma omp parallel for schedule(auto) reduction(+: count) num_threads(_num_threads)
    for (size_t i = 0; i < node_groups.size(); ++i)
    {
        const auto& node_group = node_groups[i];
        assert(node_group.size() == ghosts_per_node);
        (void)ghosts_per_node;

        /// Having a label of a node in the extended tree, we need to find the corresponding node
        /// in the original tree. We take the first ghost node, because all of them correspond to
        /// the same original node
        const auto original_node_preorder_id = _extended_mapping[node_group[0]];
        const phylo_node* original_node = *original_tree.get_by_preorder_id(original_node_preorder_id);
        const auto original_node_postorder_id = original_node->get_postorder_id();
        node_postorder_ids[i] = original_node_postorder_id;

        /// Get submatrices of probabilities for a group
        proba_group submatrices = get_submatrices(probas, node_group);
        size_t branch_count = 0;

        /// Explore k-mers of a group and store results in a hash map
        branch_hash_map group_hash_map;
        std::tie(group_hash_map, branch_count) = explore_group(submatrices);

        /// Save a hash map on disk
        save_hash_map(group_hash_map, group_hashmap_file(original_node_postorder_id));

        count += branch_count;
    }

    /// Load hash maps and merge them
    for (const auto& group_id : node_postorder_ids)
    {
        const auto hash_map = load_hash_map(group_hashmap_file(group_id));
        for (const auto& [key, score] : hash_map)
        {
            _phylo_kmer_db.insert(key, { group_id, score });
        }
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

    const auto threshold = core::score_threshold(_omega, _kmer_size);

    for (auto node_entry_ref : group)
    {
        const auto& node_entry = node_entry_ref.get();
        for (auto window = node_entry.begin(_kmer_size, threshold); window != node_entry.end(); ++window)
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

namespace boost
{
    namespace serialization
    {
        /// Serialize a hash map
        template<class Archive>
        inline void save(Archive& ar, const ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
        {
            size_t map_size = map.size();
            ar & map_size;

            for (const auto& [key, score] : map)
            {
                ar & key & score;
            }
        }

        /// Deserialize a hash map
        template<class Archive>
        inline void load(Archive& ar, ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
        {
            size_t map_size = 0;
            ar & map_size;

            for (size_t i = 0; i < map_size; ++i)
            {
                ::core::phylo_kmer::key_type key = ::core::phylo_kmer::nan_key;
                ::core::phylo_kmer::score_type score = ::core::phylo_kmer::nan_score;
                ar & key & score;
                map[key] = score;
            }
        }

        // split non-intrusive serialization function member into separate
        // non intrusive save/load member functions
        template<class Archive>
        inline void serialize(Archive& ar, ::db_builder::branch_hash_map& map, const unsigned int file_version)
        {
            boost::serialization::split_free(ar, map, file_version);
        }
    }
}

void db_builder::save_hash_map(const branch_hash_map& map, const std::string& filename) const
{
    std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa & map;
}

/// \brief Loads a hash map from file
db_builder::branch_hash_map db_builder::load_hash_map(const std::string& filename) const
{
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);

    ::db_builder::branch_hash_map map;
    ia & map;
    return map;
}


namespace rappas
{
    phylo_kmer_db build(const std::string& working_directory, const std::string& ar_probabilities_file,
                        const std::string& original_tree_file, const std::string& extended_tree_file,
                        const std::string& extended_mapping_file, const std::string& artree_mapping_file,
                        size_t kmer_size, core::phylo_kmer::score_type omega, size_t num_threads)
    {
        db_builder builder(working_directory, ar_probabilities_file, original_tree_file, extended_tree_file,
                           extended_mapping_file, artree_mapping_file, kmer_size, omega, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
