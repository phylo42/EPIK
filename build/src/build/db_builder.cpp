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
#include "pp_matrix/phyml.h"
#include "utils/utils.h"

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
                                   size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter,
                                   double mu, size_t num_threads);
    public:
        /// Member types
        /// \brief A hash map to store all the phylo-kmers, placed to one original node
        using branch_hash_map = hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;

        /// \brief A group of node ids that must be processed together. We group together
        ///        the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(const string& working_directory, const string& ar_probabilities_file,
                   const string& original_tree_file, const string& extended_tree_file,
                   const string& extended_mapping_file, const string& artree_mapping_file,
                   size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter,
                   double mu, size_t num_threads);
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
                       size_t kmer_size, xpas::phylo_kmer::score_type omega,
                       filter_type filter, double mu, size_t num_threads)
    : _working_directory{working_directory}
    ,_hashmaps_directory{(fs::path{working_directory} / fs::path{"hashmaps"}).string()}
    , _ar_probabilities_file{ar_probabilities_file}
    , _original_tree_file{original_tree_file}
    , _extended_tree_file{extended_tree_file}
    , _extended_mapping_file{extended_mapping_file}
    , _artree_mapping_file{artree_mapping_file}
    , _kmer_size{kmer_size}
    , _omega{omega}
    , _filter{filter}
    , _mu{mu}
    , _num_threads{num_threads}
    /// I do not like reading a file here, but it seems to be better than having something like "set_tree"
    /// in the public interface of the phylo_kmer_db class.
    , _phylo_kmer_db{kmer_size, omega, xpas::io::read_as_string(original_tree_file)}
{
}

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

    /// Load PhyML output
    const auto proba_matrix = rappas::io::load_phyml_probas(_ar_probabilities_file);

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

    if (_filter == filter_type::entropy)
    {
        std::cout << "Filtering by entropy" << std::endl;
    }
    else if (_filter == filter_type::max_deviation)
    {
        std::cout << "Filtering by max deviation" << std::endl;
    }
    else if (_filter == filter_type::max_difference)
    {
        std::cout << "Filtering by max difference..." << std::endl;
    }
    else if (_filter == filter_type::random)
    {
        std::cout << "Random filtering..." << std::endl;
    }
    else
    {
        std::cout << "No filtering." << std::endl;
    }

    if (_filter == filter_type::entropy || _filter == filter_type::max_deviation || _filter == filter_type::max_difference)
    {
        phylo_kmer_db temp_db(_phylo_kmer_db.kmer_size(), _phylo_kmer_db.omega(), "");

        std::cout << "Calculating entropy..." << std::endl;
        /// Load hash maps and merge them
        for (const auto group_id : group_ids)
        {
            const auto hash_map = load_hash_map(group_hashmap_file(group_id));
            for (const auto& [key, score] : hash_map)
            {
                temp_db.insert(key, {group_id, score});
            }
        }

        const auto tree = xpas::io::parse_newick(_phylo_kmer_db.tree());
        const auto threshold = xpas::score_threshold(_phylo_kmer_db.omega(), _phylo_kmer_db.kmer_size());


        hash_map<phylo_kmer::key_type, double> filter_stats;
        for (const auto& [key, entries] : temp_db)
        {
            /// calculate the score sum to normalize scores
            double score_sum = 0;
            for (const auto& [_, log_score] : entries)
            {
                score_sum += logscore_to_score(log_score);
            }

            /// do not forget the branches that are not stored in the database,
            /// they suppose to have the threshold score
            score_sum += static_cast<float>(tree.get_node_count() - entries.size()) * threshold;

            /// Entropy
            if (_filter == filter_type::entropy)
            {
                for (const auto& [branch, log_score] : entries)
                {
                    const auto weighted_score = logscore_to_score(log_score) / score_sum;
                    const auto target_value = shannon(weighted_score);
                    const auto weighted_threshold = threshold / score_sum;
                    const auto target_threshold = shannon(weighted_threshold);

                    if (filter_stats.find(key) == filter_stats.end())
                    {
                        filter_stats[key] = static_cast<float>(tree.get_node_count()) * target_threshold;
                    }
                    filter_stats[key] = filter_stats[key] - target_threshold + target_value;
                }
            }
            /// Maximum deviation
            else if (_filter == filter_type::max_deviation)
            {
                double mean_score = 0;
                for (const auto& [branch, log_score] : entries)
                {
                    double score = logscore_to_score(log_score) / score_sum;
                    mean_score += score;
                }
                mean_score /= entries.size();

                double max_diff = 0;
                for (const auto& [branch, log_score] : entries)
                {
                    double score = std::min(std::pow(10, log_score), 1.0) / score_sum;
                    max_diff = std::max(max_diff, std::abs(score - mean_score));
                }

                filter_stats[key] = max_diff;
            }
            else if (_filter == filter_type::max_difference)
            {
                double min_score = std::numeric_limits<double>::max();
                double max_score = std::numeric_limits<double>::min();
                for (const auto& [branch, log_score] : entries)
                {
                    double score = logscore_to_score(log_score) / score_sum;
                    if (score < min_score)
                    {
                        min_score = score;
                    }
                    if (score > max_score)
                    {
                        max_score = score;
                    }
                }
                filter_stats[key] = std::abs(max_score - min_score);
            }
        }

        for (const auto& [key, stats] : filter_stats)
        {
            //std::cout << stats << " ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Filtering phylo k-mers..." << std::endl;

        /// copy entropy values into a vector
        std::vector<phylo_kmer::score_type> filter_values;
        filter_values.reserve(filter_stats.size());
        for (const auto& entry : filter_stats)
        {
            filter_values.push_back(entry.second);
        }

        size_t counter = 0;
        if (_filter == filter_type::entropy)
        {
            std::cout << "Partial sort..." << std::endl;
            const auto qth_element = filter_values.size() * _mu;
            std::nth_element(filter_values.begin(), filter_values.begin() + qth_element, filter_values.end());
            const auto quantile = filter_values[qth_element];

            std::cout << "Filter threshold value: " << quantile << std::endl;

            for (const auto& [key, value] : filter_stats)
            {
                if (value <= quantile)
                {
                    ++counter;
                }
            }
            std::cout << counter << " out of " << filter_values.size() <<
                      " (" << ((float)counter) / filter_values.size() << ")" << std::endl;

            std::cout << "Merging hash maps..." << std::endl;
            for (const auto group_id : group_ids)
            {
                const auto hash_map = load_hash_map(group_hashmap_file(group_id));
                for (const auto& [key, score] : hash_map)
                {
                    if (filter_values[key] <= quantile)
                    {
                        _phylo_kmer_db.insert(key, {group_id, score});
                        //std::cout << filter_values[key] << " ";
                    }
                }
            }
            std::cout << std::endl;

        }
        else if (_filter == filter_type::max_deviation || _filter == filter_type::max_difference)
        {
            const auto qth_element = filter_values.size() * (1 - _mu);
            std::nth_element(filter_values.begin(), filter_values.begin() + qth_element, filter_values.end());
            const auto quantile = filter_values[qth_element];

            std::cout << "Filter threshold value: " << quantile << std::endl;

            for (const auto& [key, value] : filter_stats)
            {
                if (value >= quantile)
                {
                    ++counter;
                }
            }

            std::cout << counter << " out of " << filter_values.size() <<
                      " (" << ((float)counter) / filter_values.size() << ")" << std::endl;

            std::cout << "Merging hash maps..." << std::endl;
            for (const auto group_id : group_ids)
            {
                const auto hash_map = load_hash_map(group_hashmap_file(group_id));
                for (const auto& [key, score] : hash_map)
                {
                    if (filter_values[key] >= quantile)
                    {
                        _phylo_kmer_db.insert(key, {group_id, score});
                        //std::cout << score << " ";
                    }
                }
            }
            std::cout << std::endl;
        }
    }
    else if (_filter == filter_type::random)
    {
        /// create a hash map which contains boolean values
        /// which tell if we keep corresponding k-mers.
        hash_map<phylo_kmer::key_type, bool> keep_keys;
        for (const auto group_id : group_ids)
        {
            const auto hash_map = load_hash_map(group_hashmap_file(group_id));
            for (const auto& [key, _] : hash_map)
            {
                keep_keys[key] = true;
            }
        }

        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0, 1);

        /// update keep_keys randomly
        for (auto& [key, _] : keep_keys)
        {
            keep_keys[key] = distribution(generator) <= _mu;
        }

        /// Load hash maps and merge them
        std::cout << "Merging hash maps..." << std::endl;
        for (const auto group_id : group_ids)
        {
            const auto hash_map = load_hash_map(group_hashmap_file(group_id));
            for (const auto& [key, score] : hash_map)
            {
                if (keep_keys[key])
                {
                    _phylo_kmer_db.insert(key, {group_id, score});
                }
            }
        }

        size_t counter = 0;
        for (auto& [_, keep] : keep_keys)
        {
            if (keep)
            {
                counter++;
            }
        }
        std::cout << counter << " out of " << keep_keys.size() <<
                  " (" << ((float)counter) / keep_keys.size() << ")" << std::endl;
    }
    else
    {
        std::cout << "Merging hash maps..." << std::endl;
        for (const auto group_id : group_ids)
        {
            const auto hash_map = load_hash_map(group_hashmap_file(group_id));
            for (const auto& [key, score] : hash_map)
            {
                _phylo_kmer_db.insert(key, {group_id, score});
            }
        }
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
            std::cerr << "Internal error: could not find " << artree_node_label << " node." << std::endl;
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

#ifdef NDEBUG
    #pragma omp parallel for schedule(auto) reduction(+: count) default(none) num_threads(_num_threads) \
        shared(node_groups, original_tree, ghosts_per_node, node_postorder_ids, probas)
#endif
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

    const auto threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));

    for (auto node_entry_ref : group)
    {
        const auto& node_entry = node_entry_ref.get();
        //std::cout << std::endl << "NODE " << node_entry.get_label() << std::endl;
        for (auto window = node_entry.begin(_kmer_size, threshold); window != node_entry.end(); ++window)
        {
            //std::cout << std::endl << "WINDOW " << window->get_start_pos() << std::endl;
            for (const auto& kmer : *window)
            {
                //std::cout << kmer.key << " " << core::decode_kmer(kmer.key, _kmer_size) << std::endl;
                put(hash_map, kmer.key, kmer.score);
                ++count;
            }
        }
    }

    return {std::move(hash_map), count};
}

namespace boost::serialization
{
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
    inline void load(Archive& ar, ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::nan_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::nan_score;
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
                        size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter, double mu, size_t num_threads)
    {
        db_builder builder(working_directory, ar_probabilities_file, original_tree_file, extended_tree_file,
                           extended_mapping_file, artree_mapping_file, kmer_size, omega, filter, mu, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
