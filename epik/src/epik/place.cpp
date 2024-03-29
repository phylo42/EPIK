#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <iostream>
#include <i2l/seq.h>
#include <i2l/phylo_kmer_db.h>
#include <i2l/phylo_tree.h>
#include <i2l/kmer_iterator.h>
#include <i2l/seq_record.h>
#include <i2l/fasta.h>
#include <epik/place.h>

#include <chrono>

#if defined(EPIK_OMP)
#include <omp.h>
#endif

#if defined(EPIK_SSE) or \
    defined(EPIK_AVX) or \
    defined(EPIK_AVX2) or \
    defined(EPIK_AVX512)
#include <epik/intrinsic.h>
#endif


using namespace epik::impl;
using namespace epik;
using i2l::seq_record;


#ifdef __clang__

#include <cmath>

namespace epik::impl
{
    float (*pow)(float, float) = std::pow;
}

#else

namespace epik::impl
{
    double (*pow)(double, double) = std::pow;
}
#endif


/// \brief Copies the keys of an input map to a vector
std::vector<std::string_view> copy_keys(const sequence_map_t& map)
{
    std::vector<std::string_view> values;
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

placer::placer(const i2l::phylo_kmer_db& db, const i2l::phylo_tree& original_tree,
               size_t keep_at_most, double keep_factor, size_t num_threads)
    : _db{ db }
    , _original_tree{ original_tree }
    , _threshold{ i2l::score_threshold(db.omega(), db.kmer_size()) }
    , _log_threshold{ std::log10(_threshold) }
    , _keep_at_most{ keep_at_most }
    , _keep_factor{ keep_factor }
    , _max_threads{ std::max(num_threads, 1ul) }
    , _scores(_max_threads, score_vector(original_tree.get_node_count()))
    , _scores_amb(_max_threads, score_vector(original_tree.get_node_count()))
    , _counts(_max_threads, count_vector(original_tree.get_node_count()))
    , _counts_amb(_max_threads, count_vector(original_tree.get_node_count()))
    , _edges(_max_threads)
{
    /// precompute pendant lengths
    for (i2l::phylo_kmer::branch_type i = 0; i < original_tree.get_node_count(); ++i)
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

bool compare_placed_branches(const placement& lhs, const placement& rhs)
{
    return lhs.score > rhs.score;
}

/// \brief Selects keep_at_most most placed branches among these that have count > 0
std::vector<placement> placer::select_best_placements(std::vector<placement> placements, size_t num_kmers)
{
    /// Partially select best keep_at_most placements
    size_t return_size = std::min(_keep_at_most, placements.size());

    /// if no single query k-mer was found, all counts are zeros, and
    /// we create first keep_at_most placements
    if (return_size == 0)
    {
        return_size = _keep_at_most;
        placements.reserve(return_size);

        const auto threshold_score = _log_threshold * static_cast<i2l::phylo_kmer::score_type>(num_kmers)
            / static_cast<i2l::phylo_kmer::score_type>(_db.kmer_size());
        for (size_t i = 0; i < _keep_at_most; ++i)
        {
            placements.push_back({ i2l::phylo_kmer::branch_type (i), threshold_score, 0.0, 0, 0.0, 0.0 });
        }
    }
    std::partial_sort(std::begin(placements),
                      std::begin(placements) + (long)return_size,
                      std::end(placements),
                      compare_placed_branches);
    placements.resize(return_size);
    return placements;
}


/// \brief Transforms (pow10) the scores of all placements from an input array and sums it up
/// We use a longer float type, not phylo_kmer::score_type here, because 10 ** score can be a small number.
placement::weight_ratio_type placer::sum_scores(const std::vector<placement>& placements, std::string_view seq)
{
    const auto num_branches = static_cast<i2l::phylo_kmer::score_type>(_original_tree.get_node_count());
    const auto num_placements = static_cast<i2l::phylo_kmer::score_type>(placements.size());
    const auto num_kmers = static_cast<i2l::phylo_kmer::score_type>(seq.size() - _db.kmer_size() + 1);
    const auto kmer_size = static_cast<i2l::phylo_kmer::score_type>(_db.kmer_size());

    /// There are n branches where we placed the sequence, and N-n where we did not.
    /// We score each of them with (#kmers * log_threshold) / k. This gives us the total score for
    /// all branches to which the query was not placed:
    placement::weight_ratio_type sum_not_placed =
        (num_branches - num_placements) * epik::impl::pow(10.0, (num_kmers * _log_threshold / kmer_size));

    /// The final sum includes branches that were scored by k-mers explicitly
    placement::weight_ratio_type sum_placed = 0.0f;
    for (const auto& placement : placements)
    {
        sum_placed += epik::impl::pow(10.0, placement::weight_ratio_type(placement.score));
    }
    return sum_not_placed + sum_placed;
}

/// \brief Copies placements that have a weight ratio >= some threshold value. The threshold
/// is calculated as a relative _keep_factor from a maximum weight_ratio among the given placements.
std::vector<placement> filter_by_ratio(const std::vector<placement>& placements, double _keep_factor)
{
    /// calculate the ratio threshold. Here we assume that input placements are sorted
    const auto best_ratio = placements.empty() ? 0.0f : placements[0].weight_ratio;
    const auto ratio_threshold = best_ratio * _keep_factor;

    std::vector<placement> result;
    result.reserve(placements.size());
    std::copy_if(std::begin(placements), std::end(placements), std::back_inserter(result),
                 [ratio_threshold](const placement& p) { return p.weight_ratio >= ratio_threshold; });
    return result;
}

placed_collection placer::place(const std::vector<seq_record>& seq_records, size_t num_threads)
{
    (void)num_threads;

    /// There may be identical sequences with different headers. We group them
    /// by the sequence content to not to place the same sequences more than once
    const auto sequence_map = group_by_sequence_content(seq_records);

    /// To support OpenMP, we need to iterate over unique sequences in the old-style fashion.
    /// To do this, we copy all the unique keys from a map to a vector.
    /// Keys are std::string_view's, so copying is cheap enough
    const auto unique_sequences = copy_keys(sequence_map);

    /// Place only unique sequences
    std::vector<placed_sequence> placed_seqs(unique_sequences.size());

    //const auto begin_omp = std::chrono::steady_clock::now();
#ifdef EPIK_OMP
    #if __clang__
    #pragma omp parallel for schedule(dynamic) num_threads(num_threads) \
        default(none) shared(epik::impl::pow, unique_sequences, placed_seqs)
    #elif defined (__GNUC__) && (__GNUC__ < 9)
    /// In pre-GCC-9, const variables (unique_sequences here) were predefined shared automatically
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) default(none) shared(epik::impl::pow, placed_seqs)
    #else
#pragma omp parallel for schedule(dynamic) num_threads(num_threads) \
    default(none) shared(epik::impl::pow, unique_sequences, placed_seqs)
    #endif
#endif
    for (size_t i = 0; i < unique_sequences.size(); ++i)
    {
        auto keep_factor = _keep_factor;
        const auto sequence = unique_sequences[i];

        placed_seqs[i] = place_seq(sequence);

        /// compute weight ratio
        const auto score_sum = sum_scores(placed_seqs[i].placements, sequence);
        const auto num_kmers = sequence.size() - _db.kmer_size() + 1;
        placed_seqs[i].placements = select_best_placements(std::move(placed_seqs[i].placements), num_kmers);
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
                const auto power = epik::impl::pow(10.0f, placement::weight_ratio_type(placement.score));
                if (power == 0.0)
                {
                    placement.weight_ratio = 0.0;
                }
                else
                {
                    placement.weight_ratio = power / score_sum;
                }
            }
        }

        /// Remove placements with low weight ratio
        placed_seqs[i].placements = filter_by_ratio(placed_seqs[i].placements, keep_factor);
    }
    //const auto end_omp = std::chrono::steady_clock::now();
    //const float seconds = (float)std::chrono::duration_cast<std::chrono::milliseconds>(
    //    end_omp - begin_omp).count() / 1000.0f;
    //float speed = seq_records.size() / seconds / num_threads;
    //std::cout << "Query/sec (per thread): " << speed << std::endl << std::endl;
    return { sequence_map, std::move(placed_seqs) };
}


auto query_kmers(std::string_view seq, const i2l::phylo_kmer_db& db)
{
    using search_result = std::vector<decltype(db.search(0))>;

    /// Results of DB search for exact k-mers and
    /// pairs (k-mer, results of search) for ambiguous k-mers
    struct kmer_results
    {
        search_result exact;
        std::vector<search_result> ambiguous;
    };
    kmer_results result;

    result.exact.reserve(seq.size() - db.kmer_size() + 1);

    /// Query every k-mer that has no more than one ambiguous character
    for (const auto& [kmer, keys] : i2l::to_kmers<i2l::one_ambiguity_policy>(seq, db.kmer_size()))
    {
        (void) kmer;
        if (keys.size() == 1)
        {
            const auto key = keys[0];
            auto key_result = db.search(key);
            if (key_result)
            {
                result.exact.push_back(key_result);
            }
        }
        else
        {
            for (const auto& key : keys)
            {
                result.ambiguous.emplace_back();
                result.ambiguous.back().push_back(db.search(key));
            }
        }
    }
    return result;
}


/// \brief Places a fasta sequence
placed_sequence placer::place_seq(std::string_view seq)
{
    const auto num_of_kmers = seq.size() - _db.kmer_size() + 1;

#if defined(EPIK_OMP)
    const auto thread_id = omp_get_thread_num();
#else
    const size_t thread_id = 0;
#endif
    auto& thread_counts = _counts[thread_id];
    auto& thread_scores = _scores[thread_id];
    auto& thread_counts_amb = _counts_amb[thread_id];
    auto& thread_scores_amb = _scores_amb[thread_id];
    auto& thread_edges = _edges[thread_id];

    for (const auto& edge : thread_edges)
    {
        thread_counts[edge] = 0;
        thread_scores[edge] = 0.0f;
        thread_counts_amb[edge] = 0;
        thread_scores_amb[edge] = 0.0f;
    }
    thread_edges.clear();

    /// Let's query every k-mer in advance. We'll apply the scores later
    const auto search_results = query_kmers(seq, _db);
    const auto exact_phylo_kmers = search_results.exact;

    /// Now let's update the score vectors according to retrieved values
    for (auto exact_result : exact_phylo_kmers)
    {
        if (exact_result)
        {

#if defined(EPIK_SSE) or defined(EPIK_AVX) or defined(EPIK_AVX2) or defined(EPIK_AVX512)
            update_vector(thread_scores, thread_counts, thread_edges, *exact_result);
#else

            for (const auto& [postorder_node_id, score] : *exact_result)
            {
                if (thread_counts[postorder_node_id] == 0)
                {
                    thread_edges.push_back(postorder_node_id);
                }

                ++thread_counts[postorder_node_id];
                thread_scores[postorder_node_id] += score;
            }
#endif
        }

    }

    const auto ambiguous_phylo_kmers = search_results.ambiguous;
    /// Now let's update the score vectors according to retrieved values
    for (const auto& ambiguous_result : ambiguous_phylo_kmers)
    {
        /// hash set of branch ids that are scored by the ambiguous k-mer
        std::unordered_set<i2l::phylo_kmer::branch_type> l_amb;
        for (auto exact_result : ambiguous_result)
        {
            if (exact_result)
            {
                for (const auto& [postorder_node_id, score] : *exact_result)
                {
                    if (thread_counts_amb[postorder_node_id] == 0)
                    {
                        l_amb.insert(postorder_node_id);
                    }

                    thread_counts_amb[postorder_node_id] += 1;
                    thread_scores_amb[postorder_node_id] += static_cast<i2l::phylo_kmer::score_type>(std::pow(10, score));
                }

                /// Number of keys resolved from the k-mer
                const size_t w_size = _db.kmer_size();

                /// Calculate average scores
                for (const auto postorder_node_id: l_amb)
                {
                    const auto average_prob = (thread_scores_amb[postorder_node_id] +
                                               static_cast<float>(w_size - thread_counts_amb[postorder_node_id]) * _threshold)
                                              / static_cast<float>(w_size);

                    if (thread_counts[postorder_node_id] == 0)
                    {
                        thread_edges.push_back(postorder_node_id);
                    }

                    thread_counts[postorder_node_id] += 1;
                    thread_scores[postorder_node_id] += average_prob;
                }
            }
        }

    }

    /// Score correction
    for (const auto& edge: thread_edges)
    {
        thread_scores[edge] += static_cast<i2l::phylo_kmer::score_type>(num_of_kmers - thread_counts[edge]) * _log_threshold;
        thread_scores[edge] /= static_cast<i2l::phylo_kmer::score_type>(_db.kmer_size());
    }

    std::vector<placement> placements;
    placements.reserve(thread_edges.size());

    for (const auto& edge: thread_edges)
    {
        const auto node = _original_tree.get_by_postorder_id((i2l::phylo_node::id_type)edge);
        if (!node)
        {
            throw std::runtime_error("Could not find node by post-order id: " + std::to_string(edge));
        }

        const auto distal_length = (*node)->get_branch_length() / 2;
        placements.push_back({ edge, thread_scores[edge], 0.0,
                               thread_counts[edge], distal_length, _pendant_lengths[edge] });
    }
    return { seq, std::move(placements) };
}
