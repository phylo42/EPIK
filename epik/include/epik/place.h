#ifndef EPIK_PLACE_H
#define EPIK_PLACE_H

#include <vector>
#include <unordered_map>
#include <i2l/phylo_kmer.h>
#include <i2l/phylo_kmer_db.h>
#include <i2l/phylo_tree.h>

#ifdef __clang__
/// Clang still does not fully support boost::multiprecion.
/// We use floats instead despite of known precion-related
/// issues in rare cases

#include <cmath>

namespace epik::impl
{
    using lwr_type = float;
    constexpr auto pow = powf;
}
#else
#include <boost/multiprecision/float128.hpp>

namespace epik::impl
{
    using lwr_type = boost::multiprecision::float128;
    auto pow = [](const lwr_type& base, const lwr_type& exponent) {
        return boost::multiprecision::pow(base, exponent);
    };
}
#endif

namespace i2l
{
    class seq_record;
}

namespace epik::impl
{
    /// A mapping "sequence content -> list of headers" to group identical reads
    using sequence_map_t = std::unordered_map<std::string_view, std::vector<std::string_view>>;

    /// A placement of one sequence
    struct placement {
    public:
        using weight_ratio_type = lwr_type;

        i2l::phylo_kmer::branch_type branch_id;
        i2l::phylo_kmer::score_type score;
        weight_ratio_type weight_ratio;
        size_t count;
        i2l::phylo_node::branch_length_type distal_length;
        i2l::phylo_node::branch_length_type pendant_length;
    };

    /// A wrapper to store a sequence and its placement information
    struct placed_sequence {
        std::string_view sequence;
        std::vector<placement> placements;
    };

    /// \brief A collection of placed sequences
    /// \details Keys of the map must correspond to the sequences of the vector.
    struct placed_collection {
        impl::sequence_map_t sequence_map;
        std::vector<placed_sequence> placed_seqs;
    };
}

namespace epik
{
    /// \brief Places a collection of fasta sequences
    class placer
    {
        using placed_collection = impl::placed_collection;
        using placed_sequence = impl::placed_sequence;

    public:
        /// \brief Constructor.
        /// \details: WARNING: db and tree are stored as references to avoid copying and
        /// the overhead of smart pointers. Make sure that the lifetime of these variables is
        /// longer than placer's one.
        placer(const i2l::phylo_kmer_db& db, const i2l::phylo_tree& _original_tree,
               size_t keep_at_most, double keep_factor);
        placer(const placer&) = delete;
        placer(placer&&) = delete;
        placer& operator=(const placer&) = delete;
        placer& operator=(placer&&) = delete;
        ~placer() noexcept = default;

        /// \brief Places a collection of fasta sequences
        placed_collection place(const std::vector<i2l::seq_record>& seq_records, size_t num_threads);

    private:

        /// \brief Places a fasta sequence
        placed_sequence place_seq(std::string_view seq);

        const i2l::phylo_kmer_db& _db;
        const i2l::phylo_tree& _original_tree;
        const i2l::phylo_kmer::score_type _threshold;
        const i2l::phylo_kmer::score_type _log_threshold;
        const size_t _keep_at_most;
        const double _keep_factor;

        // Corresponds to S[], the array storing the score of the query for each edge x of refT
        std::vector<i2l::phylo_kmer::score_type> _scores;
        std::vector<i2l::phylo_kmer::score_type> _scores_amb;

        // Corresponds to C[], the array counting the number of k-mers in the query mapped to x
        std::vector<size_t> _counts;
        std::vector<size_t> _counts_amb;

        /// Corresponds to L[], the list of edges mapped to some k-mer in the query
        std::vector<i2l::phylo_kmer::branch_type> _edges;

        std::vector<double> _pendant_lengths;
    };
}

#endif
