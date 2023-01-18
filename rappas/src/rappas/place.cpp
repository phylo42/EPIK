#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <xcl/phylo_kmer_db.h>
#include <xcl/phylo_tree.h>
#include <xcl/kmer_iterator.h>
#include <xcl/seq_record.h>
#include <xcl/fasta.h>
#include <iostream>
#include "place.h"

using namespace rappas::impl;
using namespace rappas;
using xcl::seq_record;


/// \brief Copies the keys of an input map to a vector
template<typename K, typename V, template<class, class> typename Map>
std::vector<K> copy_keys(const Map<K, V>& map)
{
    std::vector<K> values;
    values.reserve(map.size());

    for(const auto& [key, value] : map)
    {
        (void)value;
        values.push_back(key);
    }
    return values;
}

/// \brief Groups fasta sequences by their sequence content.
/// \details Returns a map
///     {
///       sequence_1 -> [list of read ids],
///       sequence_2 -> [list of read ids],
///       ...
///     }
/// to store identical reads together
sequence_map_t group_by_sequence_content(const std::vector<seq_record>& seq_records)
{
    sequence_map_t sequence_map;
    for (const auto& seq_record : seq_records)
    {
        sequence_map[seq_record.sequence()].push_back(seq_record.header());
    }
    return sequence_map;
}

placer::placer(const xcl::phylo_kmer_db& db, const xcl::phylo_tree& original_tree, size_t keep_at_most, double keep_factor)
    : _db{ db }
    , _original_tree{ original_tree }
    , _threshold{ xcl::score_threshold(db.omega(), db.kmer_size()) }
    , _log_threshold{ std::log10(_threshold) }
    , _keep_at_most{ keep_at_most }
    , _keep_factor{ keep_factor }
    , _scores(original_tree.get_node_count())
    , _scores_amb(original_tree.get_node_count())
    , _counts(original_tree.get_node_count())
    , _counts_amb(original_tree.get_node_count())
{
    /// precompute pendant lengths
    for (xcl::phylo_kmer::branch_type i = 0; i < original_tree.get_node_count(); ++i)
    {
        /// i is a post-order node id here. The phylo_kmer_db::search returns the post-order ids,
        /// not the pre-order ones
        const auto node = _original_tree.get_by_postorder_id(i);
        if (!node)
        {
            const auto node_id_str = std::to_string(i);
            throw std::runtime_error("Could not find node by post-order id: " + node_id_str);
        }

        const auto distal_length = (*node)->get_branch_length() / 2;

        /// For pendant_length
        const auto num_subtree_nodes = _db.tree_index()[i].subtree_num_nodes;
        const auto subtree_branch_length = _db.tree_index()[i].subtree_total_length;

        /// calculate the mean branch length in the subtree (excluding the branch with this post-order id)
        auto mean_subtree_branch_length = 0.0;
        if (num_subtree_nodes > 1)
        {
            mean_subtree_branch_length = subtree_branch_length / num_subtree_nodes;
        }

        const auto pendant_length = mean_subtree_branch_length + distal_length;
        _pendant_lengths.push_back(pendant_length);
    }
}

/// \brief Transforms (pow10) the scores of all placements from an input array and sums it up
/// We use a longer float type, not phylo_kmer::score_type here, because 10 ** score can be a small number.
placement::weight_ratio_type sum_scores(const std::vector<placement>& placements)
{
    placement::weight_ratio_type sum = 0.0;
    for (const auto& placement : placements)
    {
        sum += boost::multiprecision::pow(10.0, placement::weight_ratio_type(placement.score));
    }
    return sum;
}

/// \brief Copies placements that have a weight ratio >= some threshold value. The threshold
/// is calculated as a relative _keep_factor from a maximum weight_ratio among the given placements.
std::vector<placement> filter_by_ratio(const std::vector<placement>& placements, double _keep_factor)
{
    /// calculate the ratio threshold. Here we assume that input placements are sorted
    const auto best_ratio = placements.size() > 0 ? placements[0].weight_ratio : 0.0f;
    const auto ratio_threshold = best_ratio *_keep_factor;

    std::vector<placement> result;
    result.reserve(placements.size());
    std::copy_if(std::begin(placements), std::end(placements), std::back_inserter(result),
                 [ratio_threshold](const placement& p) { return p.weight_ratio >= ratio_threshold; });
    return result;
}

placed_collection placer::place(const std::vector<seq_record>& seq_records, size_t num_threads)
{
    /// There may be identical sequences with different headers. We group them
    /// by the sequence content to not to place the same sequences more than once
    const auto sequence_map = group_by_sequence_content(seq_records);

    /// To support OpenMP, we need to iterate over unique sequences in the old-style fashion.
    /// To do this, we copy all the unique keys from a map to a vector.
    /// Keys are std::string_view's, so copying is cheap enough
    const auto unique_sequences = copy_keys(sequence_map);


    /// Place only unique sequences
    std::vector<placed_sequence> placed_seqs(unique_sequences.size());
    (void)num_threads;
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads) \
        default(none) shared(placed_seqs, unique_sequences)
    for (size_t i = 0; i < unique_sequences.size(); ++i)
    {
        auto keep_factor = _keep_factor;
        const auto sequence = unique_sequences[i];

        placed_seqs[i] = place_seq(sequence);

        /// compute weight ratio
        const auto score_sum = sum_scores(placed_seqs[i].placements);
        for (auto& placement : placed_seqs[i].placements)
        {
            /// If the scores are that small that taking 10 to these powers is still zero
            /// according to boost::multiprecision::pow, then score_sum is zero.
            /// Assign all weight_ration to zeros and keep_factor to zero as well to
            /// not filter them out.
            if (score_sum == 0)
            {
                placement.weight_ratio = 0.0f;
                keep_factor = 0.0f;
            }
            else
            {
                placement.weight_ratio = boost::multiprecision::pow(10.0f, placement::weight_ratio_type(placement.score)) / score_sum;
            }
        }

        /// Remove placements with low weight ratio
        placed_seqs[i].placements = filter_by_ratio(placed_seqs[i].placements, keep_factor);
    }
    return { sequence_map, placed_seqs };
}

bool compare_placed_branches(const placement& lhs, const placement& rhs)
{
    return lhs.score > rhs.score;
}

/// \brief Selects keep_at_most most placed branches among these that have count > 0
std::vector<placement> select_best_placements(std::vector<placement> placements, size_t keep_at_most)
{
    /// Partially select best keep_at_most placements
    size_t return_size = std::min(keep_at_most, placements.size());

    /// if no single query k-mer was found, all counts are zeros, and
    /// we take just first keep_at_most branches
    if (return_size == 0)
    {
        return_size = std::min(keep_at_most, placements.size());
        placements.resize(return_size);
        std::copy(placements.begin(), placements.end(), placements.begin());
    }
    std::partial_sort(std::begin(placements),
                      std::begin(placements) + return_size,
                      std::end(placements),
                      compare_placed_branches);
    placements.resize(return_size);
    return placements;
}

/// \brief Places a fasta sequence
placed_sequence placer::place_seq(std::string_view seq)
{
    const auto num_of_kmers = seq.size() - _db.kmer_size() + 1;

    for (const auto& edge : _edges)
    {
        _counts[edge] = 0;
        _scores[edge] = 0.0f;
        _counts_amb[edge] = 0;
        _scores_amb[edge] = 0.0f;
    }
    _edges.clear();

    /// Query every k-mer that has no more than one ambiguous character
    for (const auto& [kmer, keys] : xcl::to_kmers<xcl::one_ambiguity_policy>(seq, _db.kmer_size()))
    {
        (void)kmer;

        /// if the k-mer is unambiguous
        if (keys.size() == 1)
        {
            const auto key = keys[0];

            /// Update placements if found
            if (auto entries = _db.search(key); entries)
            {
#ifdef KEEP_POSITIONS
                for (const auto& [postorder_node_id, score, position] : *entries)
                {
                    (void)position;
#else
                for (const auto& [postorder_node_id, score] : *entries)
                {
#endif
                    if (_counts[postorder_node_id] == 0)
                    {
                        _edges.push_back(postorder_node_id);
                    }

                    _counts[postorder_node_id] += 1;
                    _scores[postorder_node_id] += score;
                }

            }
        }
        /// treat ambiguities with mean
        else
        {
            /// hash set of branch ids that are scored by the ambiguous k-mer
            std::unordered_set<xcl::phylo_kmer::branch_type> l_amb;

            for (const auto& key : keys)
            {
                if (auto entries = _db.search(key); entries)
                {
#ifdef KEEP_POSITIONS
                    for (const auto& [postorder_node_id, score, position] : *entries)
                    {
                        (void)position;
#else
                    for (const auto& [postorder_node_id, score] : *entries)
                    {
#endif
                        if (_counts_amb[postorder_node_id] == 0)
                        {
                            l_amb.insert(postorder_node_id);
                        }

                        _counts_amb[postorder_node_id] += 1;
                        _scores_amb[postorder_node_id] += static_cast<xcl::phylo_kmer::score_type>(std::pow(10, score));
                    }
                }
            }

            /// Number of keys resolved from the k-mer
            const size_t w_size = keys.size();

            /// Calculate average scores
            for (const auto postorder_node_id : l_amb)
            {
                const auto average_prob = (
                    _scores_amb[postorder_node_id] +
                    static_cast<float>(w_size - _counts_amb[postorder_node_id]) * _threshold
                    ) / static_cast<float>(w_size);

                if (_counts[postorder_node_id] == 0)
                {
                    _edges.push_back(postorder_node_id);
                }

                _counts[postorder_node_id] += 1;
                _scores[postorder_node_id] += average_prob;
            }
        }
    }

    /// Score correction
    for (const auto& edge: _edges)
    {
        _scores[edge] += static_cast<xcl::phylo_kmer::score_type>(num_of_kmers - _counts[edge]) * _log_threshold;
    }

    std::vector<placement> placements;
    for (const auto& edge: _edges)
    {
        const auto node = _original_tree.get_by_postorder_id(edge);
        if (!node)
        {
            throw std::runtime_error("Could not find node by post-order id: " + std::to_string(edge));
        }

        const auto distal_length = (*node)->get_branch_length() / 2;
        placements.push_back({edge, _scores[edge], 0.0, _counts[edge], distal_length, _pendant_lengths[edge] });
    }

    return { seq, select_best_placements(std::move(placements), _keep_at_most) };
}
