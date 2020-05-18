#include <vector>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <xpas/phylo_kmer_db.h>
#include <xpas/phylo_tree.h>
#include <xpas/kmer_iterator.h>
#include <utils/io/fasta.h>
#include <unordered_set>
#include "place.h"

using xpas::io::fasta;
using namespace rappas::impl;
using namespace rappas;


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
sequence_map_t group_by_sequence_content(const std::vector<fasta>& fasta_collection)
{
    sequence_map_t sequence_map;
    for (const auto& seq_record : fasta_collection)
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

placed_collection placer::place(const std::vector<fasta>& fasta_collection, size_t num_threads) const
{
    /// There may be identical sequences with different headers. We group them
    /// by the sequence content to not to place the same sequences more than once
    const auto sequence_map = group_by_sequence_content(fasta_collection);

    /// To support OpenMP, we need to iterate over unique sequences in the old-style fashion.
    /// To do this, we copy all the unique keys from a map to a vector.
    /// Keys are std::string_view's, so copying is cheap enough
    const auto unique_sequences = copy_keys(sequence_map);

    /// Place only unique sequences
    std::vector<placed_sequence> placed_seqs(unique_sequences.size());
    (void)num_threads;
    #pragma omp parallel for schedule(auto) num_threads(num_threads)
    for (size_t i = 0; i < unique_sequences.size(); ++i)
    {
        const auto sequence = unique_sequences[i];
        const auto headers = sequence_map.at(sequence);

        placed_seqs[i] = place_seq(sequence);

        /// compute weight ratio
        const auto score_sum = sum_scores(placed_seqs[i].placements);
        for (auto& placement : placed_seqs[i].placements)
        {
            placement.weight_ratio = boost::multiprecision::pow(10.0f, placement::weight_ratio_type(placement.score)) / score_sum;
        }

        /// Remove placements with low weight ratio
        placed_seqs[i].placements = filter_by_ratio(placed_seqs[i].placements, _keep_factor);
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
    const size_t return_size = std::min(keep_at_most, result.size());
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
                for (const auto& [postorder_node_id, score] : *entries)
                {
                    placements[postorder_node_id].count += 1;
                    placements[postorder_node_id].score += score - _log_threshold;
                }
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
                    for (const auto& [postorder_node_id, score] : *entries)
                    {
                        l_amb.insert(postorder_node_id);

                        ambiguous_placements[postorder_node_id].count += 1;
                        ambiguous_placements[postorder_node_id].score += std::pow(10, score);
                    }
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

    return { seq, select_best_placements(placements, _keep_at_most) };
}
