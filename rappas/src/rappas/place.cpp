#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <iostream>
#include <xpas/phylo_kmer_db.h>
#include <xpas/phylo_tree.h>
#include <xpas/kmer_iterator.h>
#include <xpas/seq_record.h>
#include <xpas/fasta.h>
#include "place.h"

using namespace rappas::impl;
using namespace rappas;
using xpas::seq_record;


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

placer::placer(const xpas::phylo_kmer_db& db, const xpas::phylo_tree& original_tree, size_t keep_at_most, double keep_factor) noexcept
    : _db{ db }
    , _original_tree{ original_tree }
    , _threshold{ xpas::score_threshold(db.omega(), db.kmer_size()) }
    , _log_threshold{ std::log10(_threshold) }
    , _keep_at_most{ keep_at_most }
    , _keep_factor{ keep_factor }
{}

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

placed_collection placer::place(const std::vector<seq_record>& seq_records, size_t num_threads) const
{
    /// There may be identical sequences with different headers. We group them
    /// by the sequence content to not to place the same sequences more than once
    const auto sequence_map = group_by_sequence_content(seq_records);

    /// To support OpenMP, we need to iterate over unique sequences in the old-style fashion.
    /// To do this, we copy all the unique keys from a map to a vector.
    /// Keys are std::string_view's, so copying is cheap enough
    const auto unique_sequences = copy_keys(sequence_map);

    auto keep_factor = _keep_factor;

    /// Place only unique sequences
    std::vector<placed_sequence> placed_seqs(unique_sequences.size());
    (void)num_threads;
    /*#pragma omp parallel for schedule(auto) num_threads(num_threads)*/
    for (size_t i = 0; i < unique_sequences.size(); ++i)
    {
        /// place the sequence
        placed_seqs[i] = std::move(place_seq(unique_sequences[i]));

        /// filter placements by likelihood weight ratio
        placed_seqs[i].placements = filter_by_ratio(placed_seqs[i].placements, keep_factor);
    }
    return { sequence_map, placed_seqs };
}

bool compare_placed_branches(const placement& lhs, const placement& rhs)
{
    return lhs.score > rhs.score;
}

/// \brief Selects keep_at_most most placed branches among these that have count > 0
std::vector<placement> select_best_placements(const std::vector<placement>& placements, size_t keep_at_most)
{
    /// WARNING: Here we copy the placements to make sorting easier. It may be not
    /// the most effective and elegant solution.
    std::vector<placement> result(placements.size());
    auto it = std::copy_if(placements.begin(), placements.end(), result.begin(),
                           [](const placement& pb){ return pb.count > 0; } );
    result.resize(std::distance(result.begin(), it));

    /// Partially select best keep_at_most placements
    size_t return_size = std::min(keep_at_most, result.size());

    /// if no single query k-mer was found, all counts are zeros, and
    /// we take just first keep_at_most branches
    if (return_size == 0)
    {
        return_size = std::min(keep_at_most, placements.size());
        result.resize(return_size);
        std::copy(placements.begin(), placements.end(), result.begin());

    }
    std::partial_sort(std::begin(result), std::begin(result) + return_size, std::end(result), compare_placed_branches);
    return { std::begin(result), std::begin(result) + return_size };
}

/// \brief Places a fasta sequence
placed_sequence placer::place_seq(std::string_view seq) const
{
    const auto num_of_kmers = seq.size() - _db.kmer_size() + 1;
    const auto sequence_log_threshold = num_of_kmers * _log_threshold;
    const auto num_branch_nodes = _original_tree.get_node_count();

    /// Initialize "ambiguous placements" array, which correspond to S_amb[] and C_amb[]
    /// We initialize it here and clean it at every iteration over an umbiguous k-mer
    /// for efficiency reasons.
    std::vector<placement> ambiguous_placements(num_branch_nodes);

    /// Initialize the "placements" array. This is an array of size of the number of branches
    /// in a tree, which contains triplets [branch_id, score, count].
    /// placements.score and placements.count correspond to S[] and C[] arrays from the original paper.
    std::vector<placement> placements;
    placements.reserve(num_branch_nodes);
    for (xpas::phylo_kmer::branch_type i = 0; i < num_branch_nodes; ++i)
    {
        /// i is a post-order node id here. The phylo_kmer_db::search returns the post-order ids,
        /// not the pre-order ones
        const auto node = _original_tree.get_by_postorder_id(i);
        if (!node)
        {
            throw std::runtime_error("Could not find node by post-order id: " + std::to_string(i));
        }

        const auto distal_length = (*node)->get_branch_length() / 2;

        /// For pendant_length, calculate the total branch length in the whole subtree
        xpas::phylo_node::branch_length_type total_subtree_branch_lenght = 0;
        size_t num_subtree_nodes = 0;
        for (const auto& subtree_node : xpas::visit_subtree(*node))
        {
            total_subtree_branch_lenght += subtree_node.get_branch_length();
            ++num_subtree_nodes;
        }

        /// calculate the mean branch length in the subtree (excluding the branch with this post-order id)
        auto mean_subtree_branch_length = 0.0f;
        if (num_subtree_nodes > 1)
        {
            mean_subtree_branch_length = (total_subtree_branch_lenght - (*node)->get_branch_length()) / (num_subtree_nodes - 1.0f);
        }

        const auto pendant_length = mean_subtree_branch_length + distal_length;
        placements.push_back({ i, sequence_log_threshold, 0.0, 0, distal_length, pendant_length});
    }

    /// Query every k-mer that has no more than one ambiguous character
    for (const auto& [kmer, keys] : xpas::to_kmers<xpas::one_ambiguity_policy>(seq, _db.kmer_size()))
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
                    placements[postorder_node_id].count += 1;
                    placements[postorder_node_id].score += score - _log_threshold;
                }
#else
                //std::cout << key << " " << kmer << std::endl;
                for (const auto& [postorder_node_id, score] : *entries)
                {
                    placements[postorder_node_id].count += 1;
                    placements[postorder_node_id].score += score - _log_threshold;

                    //std::cout << "\t" << postorder_node_id << " " << score << " " << std::pow(10, score) << std::endl;

                }
#endif
            }
        }
        /// treat ambiguities with mean
        else
        {
            /// hash set of branch ids that are scored by the ambiguous k-mer
            std::unordered_set<xpas::phylo_kmer::branch_type> l_amb;

            for (const auto& key : keys)
            {
                if (auto entries = _db.search(key); entries)
                {
#ifdef KEEP_POSITIONS
                    for (const auto& [postorder_node_id, score, position] : *entries)
                    {
                        (void)position;
                        l_amb.insert(postorder_node_id);

                        ambiguous_placements[postorder_node_id].count += 1;
                        ambiguous_placements[postorder_node_id].score += std::pow(10, score);
                    }
#else
                    for (const auto& [postorder_node_id, score] : *entries)
                    {
                        l_amb.insert(postorder_node_id);

                        ambiguous_placements[postorder_node_id].count += 1;
                        ambiguous_placements[postorder_node_id].score += std::pow(10, score);
                    }
#endif
                }
            }

            /// Number of keys resolved from the k-mer
            const size_t w_size = keys.size();

            /// Calculate average scores
            for (const auto postorder_node_id : l_amb)
            {
                const auto average_prob = (
                    ambiguous_placements[postorder_node_id].score +
                    (w_size - ambiguous_placements[postorder_node_id].count) * _threshold) / w_size;

                placements[postorder_node_id].count += 1;
                placements[postorder_node_id].score += std::log10(average_prob) - _log_threshold;

                /// clean the ambiguous placements array
                ambiguous_placements[postorder_node_id].count = 0;
                ambiguous_placements[postorder_node_id].score = 0;
            }
        }
    }

    /// get the score of the best placement
    const auto score_compare = [](const auto& p1, const auto& p2) { return p1.score < p2.score; };
    const auto max_score = std::max_element(placements.begin(), placements.end(),
                                            score_compare)->score;

    /// compute the sum of 10^(score - max_score)
    placement::weight_ratio_type score_sum = 0.0;
    for (const auto& placement : placements)
    {
        score_sum += std::pow(10.0, placement::weight_ratio_type(placement.score - max_score));
    }

    /// compute likelihood weight ratios
    for (auto& placement : placements)
    {
        placement.weight_ratio = std::pow(10.0f, placement::weight_ratio_type(placement.score - max_score)) / score_sum;
    }

    return { seq, select_best_placements(placements, _keep_at_most) };
}
