#include <iostream>
#include <unordered_set>
#include <chrono>
#include <random>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/phylo_tree.h>
#include <utils/io/file_io.h>
#include <xpas/newick.h>
#include "db_builder.h"
#include "pp_matrix/proba_matrix.h"
#include "pp_matrix/ar.h"

using std::string;
using std::cout, std::endl;
using std::to_string;
using namespace xpas;
namespace fs = boost::filesystem;


namespace rappas
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(const string& working_directory, const string& ar_probabilities_file,
                                   const string& original_tree_file, const string& extended_tree_file,
                                   const string& extended_mapping_file, const string& artree_mapping_file,
                                   bool merge_branches, size_t kmer_size, xpas::phylo_kmer::score_type omega,
                                   filter_type filter, double mu, size_t num_threads);
    public:
        /// Member types
        /// \brief A hash map to store all the phylo-kmers, placed to one original node
#ifdef KEEP_POSITIONS

        struct score_pos_pair
        {
            phylo_kmer::score_type score;
            phylo_kmer::pos_type position;
        };

        using branch_hash_map = hash_map<phylo_kmer::key_type, score_pos_pair>;
#else
        using branch_hash_map = hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;
#endif

        /// \brief A group of node ids that must be processed together. We group together
        ///        the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(const string& working_directory, const string& ar_probabilities_file,
                   const string& original_tree_file, const string& extended_tree_file,
                   const string& extended_mapping_file, const string& artree_mapping_file,
                   bool merge_branches, size_t kmer_size, xpas::phylo_kmer::score_type omega,
                   filter_type filter, double mu, size_t num_threads);
        db_builder(const db_builder&) = delete;
        db_builder(db_builder&&) = delete;
        db_builder& operator=(const db_builder&) = delete;
        db_builder& operator=(db_builder&&) = delete;
        ~db_builder() noexcept = default;


        /// \brief Runs the database construction
        void run();

    private:
        /// \brief Results of the first stage of the algorithm:
        using explore_groups_result = std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long>;

        /// \brief The first stage of the construction algorithm. Creates a hashmap of phylo-kmers
        ///        for every node group.
        /// \return group ids (correspond to the post-order ids in the tree),
        ///         number of tuples explored, elapsed time
        std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> construct_group_hashmaps();

        /// \brief The second stage of the construction algorithm. Combines group hashmaps
        /// \return Elapsed time
        unsigned long merge_hashmaps(const std::vector<phylo_kmer::branch_type>& group_ids);

        /// \brief Returns a filename for a hashmap of a given group
        std::string group_hashmap_file(const branch_type& group) const;

        /// \brief Groups ghost nodes by corresponding original node id
        std::vector<id_group> group_ghost_ids(const std::vector<std::string>& ghost_ids) const;

        /// \brief Groups references to submatrices of probabilites, corresponding to a group of nodes
        proba_group get_submatrices(const proba_matrix& probas, const id_group& group) const;

        /// \brief Runs a phylo-kmer exploration for every ghost node of the extended_tree
        /// \return 1) A vector of group ids, which correspond to post-order node ids in the tree
        ///         2) The number of explored phylo-kmers. This number can be more than a size of a resulting database
        std::tuple<std::vector<phylo_kmer::branch_type>, size_t> explore_kmers(const phylo_tree& original_tree,
                                                                               const phylo_tree& extended_tree,
                                                                               const proba_matrix& probas);

        /// \brief Explores phylo-kmers of a collection of ghost nodes. Here we assume that the nodes
        ///        in the group correspond to one original node
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

        bool _merge_branches;

        size_t _kmer_size;
        xpas::phylo_kmer::score_type _omega;

        filter_type _filter;
        double _mu;

        size_t _num_threads;
        phylo_kmer_db _phylo_kmer_db;

        extended_mapping _extended_mapping;
        artree_label_mapping _artree_mapping;
    };
}

using namespace rappas;

db_builder::db_builder(const string& working_directory, const string& ar_probabilities_file,
                       const string& original_tree_file, const string& extended_tree_file,
                       const string& extended_mapping_file, const string& artree_mapping_file,
                       bool merge_branches, size_t kmer_size, xpas::phylo_kmer::score_type omega,
                       filter_type filter, double mu, size_t num_threads)
    : _working_directory{ working_directory }
    ,_hashmaps_directory{ (fs::path{working_directory} / fs::path{"hashmaps"}).string() }
    , _ar_probabilities_file{ ar_probabilities_file }
    , _original_tree_file{ original_tree_file }
    , _extended_tree_file{ extended_tree_file }
    , _extended_mapping_file{ extended_mapping_file }
    , _artree_mapping_file{ artree_mapping_file }
    , _merge_branches{ merge_branches }
    , _kmer_size{ kmer_size }
    , _omega{ omega }
    , _filter{ filter }
    , _mu{ mu }
    , _num_threads{ num_threads }
    /// I do not like reading a file here, but it seems to be better than having something like "set_tree"
    /// in the public interface of the phylo_kmer_db class.
    , _phylo_kmer_db{ kmer_size, omega, xpas::io::read_as_string(original_tree_file) }
{
}

void create_directory(const std::string& dirname)
{
    /// Create if does not exist
    if (!fs::is_directory(dirname) || !fs::exists(dirname))
    {
        if (!fs::create_directory(dirname))
        {
            throw std::runtime_error("Cannot create directory " + dirname);
        }
    }
}

void db_builder::run()
{
    /// The first stage of the algorithm -- create a hashmap for every group node
    const auto& [group_ids, num_tuples, construction_time] = construct_group_hashmaps();

    /// The second stage of the algorithm -- combine hashmaps
    const auto merge_time = merge_hashmaps(group_ids);

    /// Calculate the number of phylo-kmers stored in the database
    size_t total_entries = 0;
    for (const auto& kmer_entry : _phylo_kmer_db)
    {
        total_entries += kmer_entry.second.size();
    }

    std::cout << "Built " << total_entries << " phylo-kmers out of " << num_tuples << " for "
              << _phylo_kmer_db.size() << " k-mer values.\nTime (ms): "
              << construction_time + merge_time << "\n\n" << std::flush;
}

std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> db_builder::construct_group_hashmaps()
{
    /// Load .tsv files
    _extended_mapping = rappas::io::load_extended_mapping(_extended_mapping_file);
    _artree_mapping = rappas::io::load_artree_mapping(_artree_mapping_file);

    /// create a temporary directory for hashmaps
    create_directory(_hashmaps_directory);

    /// Load .newick files
    const auto original_tree = xpas::io::load_newick(_original_tree_file);
    const auto extended_tree = xpas::io::load_newick(_extended_tree_file);

    /// Load ancestral reconstruction output
    const auto proba_matrix = rappas::io::load_ar(_ar_probabilities_file);

    /// Run the branch and bound algorithm
    std::cout << "Building database..." << std::endl;

    const auto begin = std::chrono::steady_clock::now();
    const auto& [group_ids, num_tuples] = explore_kmers(original_tree, extended_tree, proba_matrix);
    const auto end = std::chrono::steady_clock::now();
    const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    return { group_ids, num_tuples, elapsed_time };
}

double shannon(double x)
{
    return - x * std::log2(x);
}

xpas::phylo_kmer::score_type logscore_to_score(xpas::phylo_kmer::score_type log_score)
{
    return std::min(std::pow(10, log_score), 1.0);
}

unsigned long db_builder::merge_hashmaps(const std::vector<phylo_kmer::branch_type>& group_ids)
{
    const auto begin = std::chrono::steady_clock::now();

    std::cout << "Merging hash maps..." << std::endl;
    for (const auto group_id : group_ids)
    {
        const auto hash_map = load_hash_map(group_hashmap_file(group_id));

#ifdef KEEP_POSITIONS
        if (_merge_branches)
        {
            for (const auto& [key, score_pos_pair] : hash_map)
            {
                const auto& [score, position] = score_pos_pair;
                if (auto entries = _phylo_kmer_db.search(key); entries)
                {
                    /// If there are entries, there must be only one because
                    /// we always take maximum score among different branches.
                    /// So this loop will have only one iteration
                    for (const auto& [old_node, old_score, old_position] : *entries)
                    {
                        if (old_score < score)
                        {
                            _phylo_kmer_db.replace(key, { group_id, score, position });
                        }
                    }
                }
                else
                {
                    _phylo_kmer_db.unsafe_insert(key, { group_id, score, position });
                }
            }
        }
        else
        {
            for (const auto& [key, score_pos_pair] : hash_map)
            {
                const auto& [score, position] = score_pos_pair;
                _phylo_kmer_db.unsafe_insert(key, { group_id, score, position });
            }
        }
#else
        if (_merge_branches)
        {
            throw std::runtime_error("--merge-branches is only supported for xpas compiled with the KEEP_POSITIONS flag.");
        }
        else
        {
            for (const auto& [key, score] : hash_map)
            {
                _phylo_kmer_db.unsafe_insert(key, { group_id, score });
            }
        }
#endif
    }

    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
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
    return (fs::path{_hashmaps_directory} / fs::path{std::to_string(group)}).string();
}

std::vector<db_builder::id_group> db_builder::group_ghost_ids(const std::vector<std::string>& ghost_ids) const
{
    std::vector<id_group> groups;
    groups.reserve(ghost_ids.size() / 2);

    std::unordered_map<branch_type, size_t> mapping;
    for (const auto& ghost_id : ghost_ids)
    {
        const auto original_preorder_id = _extended_mapping.at(ghost_id);
        if (const auto it = mapping.find(original_preorder_id); it != mapping.end())
        {
            groups[it->second].push_back(ghost_id);
        }
        else
        {
            groups.push_back({ghost_id});
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
            throw std::runtime_error("Internal error: could not find " + artree_node_label + " node. "
                                     "Make sure it is in the ARTree_id_mapping file.");
        }
    }

    return submatrices;
}

std::tuple<std::vector<phylo_kmer::branch_type>, size_t> db_builder::explore_kmers(
    const phylo_tree& original_tree, const phylo_tree& extended_tree, const proba_matrix& probas)
{
    size_t count = 0;

    /// Here we assume that every original node corresponds to two ghost nodes.
    const size_t ghosts_per_node = 2;

    /// Filter and group ghost nodes
    const auto node_groups = group_ghost_ids(get_ghost_ids(extended_tree));

    /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
    /// in a hash map for every group separately on disk.
    std::vector<phylo_kmer::branch_type> node_postorder_ids(node_groups.size());

//#ifdef NDEBUG
//    #pragma omp parallel for schedule(auto) reduction(+: count) num_threads(_num_threads) \
//        shared(original_tree, node_postorder_ids, probas)
//#endif
    for (size_t i = 0; i < node_groups.size(); ++i)
    {
        const auto& node_group = node_groups[i];
        assert(node_group.size() == ghosts_per_node);
        (void) ghosts_per_node;

        /// Having a label of a node in the extended tree, we need to find the corresponding node
        /// in the original tree. We take the first ghost node, because all of them correspond to
        /// the same original node
        const auto original_node_preorder_id = _extended_mapping[node_group[0]];
        const phylo_node* original_node = *original_tree.get_by_preorder_id(original_node_preorder_id);
        const auto original_node_postorder_id = original_node->get_postorder_id();
        node_postorder_ids[i] = original_node_postorder_id;

        /// Get submatrices of probabilities for a group
        proba_group submatrices = get_submatrices(probas, node_group);

        /// Explore k-mers of a group and store results in a hash map
        const auto& [group_hash_map, branch_count] = explore_group(submatrices);

        /// Save a hash map on disk
        save_hash_map(group_hash_map, group_hashmap_file(original_node_postorder_id));

        count += branch_count;
    }
    return { node_postorder_ids, count };
}

/// \brief Puts a key-value pair in a hash map. Used to process branches in parallel
#ifdef KEEP_POSITIONS
void put(db_builder::branch_hash_map& map, const phylo_kmer& kmer)
{
    if (auto it = map.find(kmer.key); it != map.end())
    {
        if (it->second.score < kmer.score)
        {
            map[kmer.key] = { kmer.score, kmer.position };
        }
    }
    else
    {
        map[kmer.key] = { kmer.score, kmer.position };
    }
}

std::pair<db_builder::branch_hash_map, size_t> db_builder::explore_group(const proba_group& group)
{
    branch_hash_map hash_map;
    size_t count = 0;

    const auto threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));

    for (auto node_entry_ref : group)
    {
        const auto& node_entry = node_entry_ref.get();
        for (auto window = node_entry.begin(_kmer_size, threshold); window != node_entry.end(); ++window)
        {
            const auto position = window->get_start_pos();
            for (const auto& kmer : *window)
            {
                phylo_kmer positioned_kmer = { kmer.key, kmer.score, position };
                put(hash_map, positioned_kmer);
                ++count;
            }
        }
    }

    return {std::move(hash_map), count};
}
#else

void put(db_builder::branch_hash_map& map, const phylo_kmer& kmer)
{
    if (auto it = map.find(kmer.key); it != map.end())
    {
        if (it->second < kmer.score)
        {
            map[kmer.key] = kmer.score;
        }
    }
    else
    {
        map[kmer.key] = kmer.score;
    }
}

std::pair<db_builder::branch_hash_map, size_t> db_builder::explore_group(const proba_group& group)
{
    branch_hash_map hash_map;
    size_t count = 0;

    const auto threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));

    for (auto node_entry_ref : group)
    {
        const auto& node_entry = node_entry_ref.get();
        for (auto window = node_entry.begin(_kmer_size, threshold); window != node_entry.end(); ++window)
        {
            for (const auto& kmer : *window)
            {
                put(hash_map, kmer);
                ++count;
            }
        }
    }

    return {std::move(hash_map), count};
}

#endif

namespace boost::serialization
{

#ifdef KEEP_POSITIONS
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto& [key, score_pos_pair] : map)
        {
            const auto& [score, position] = score_pos_pair;
            ar & key & score & position;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            xpas::phylo_kmer::pos_type position = xpas::phylo_kmer::na_pos;
            ar & key & score & position;
            map[key] = { score, position };
        }
    }
#else
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto&[key, score] : map)
        {
            ar & key & score;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            ar & key & score;
            map[key] = score;
        }
    }
#endif

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int file_version)
    {
        boost::serialization::split_free(ar, map, file_version);
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
                        bool merge_branches,
                        size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter, double mu, size_t num_threads)
    {
        db_builder builder(working_directory, ar_probabilities_file, original_tree_file, extended_tree_file,
                           extended_mapping_file, artree_mapping_file, merge_branches,
                           kmer_size, omega, filter, mu, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
