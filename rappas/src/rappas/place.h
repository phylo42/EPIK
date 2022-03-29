#ifndef RAPPAS_PLACE_PLACE_H
#define RAPPAS_PLACE_PLACE_H

#include <vector>
#include <unordered_map>
#include <xpas/phylo_kmer.h>
#include <xpas/phylo_kmer_db.h>
#include <xpas/phylo_tree.h>
#include <boost/multiprecision/float128.hpp>

namespace xpas
{
    class seq_record;
}

namespace rappas::impl
{
    /// A mapping "sequence content -> list of headers" to group identical reads
    /// TODO: check if std::unordered_map is efficient enough
    using sequence_map_t = std::unordered_map<std::string_view, std::vector<std::string_view>>;

    /// A placement of one sequence
    struct placement {
    public:
        using weight_ratio_type = boost::multiprecision::float128;

        xpas::phylo_kmer::branch_type branch_id;
        xpas::phylo_kmer::score_type score;
        weight_ratio_type weight_ratio;
        size_t count;
        xpas::phylo_node::branch_length_type distal_length;
        xpas::phylo_node::branch_length_type pendant_length;
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

namespace rappas
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
        placer(const xpas::phylo_kmer_db& db, const xpas::phylo_tree& _original_tree,
               size_t keep_at_most, double keep_factor, bool exists) noexcept;
        placer(const placer&) = delete;
        placer(placer&&) = delete;
        placer& operator=(const placer&) = delete;
        placer& operator=(placer&&) = delete;
        ~placer() noexcept = default;

        /// \brief Places a collection of fasta sequences
        placed_collection place(const std::vector<xpas::seq_record>& seq_records, size_t num_threads) const;

    private:

        /// \brief Places a fasta sequence
        placed_sequence place_seq(std::string_view seq) const;

        const xpas::phylo_kmer_db& _db;
        const xpas::phylo_tree& _original_tree;
        const xpas::phylo_kmer::score_type _threshold;
        const xpas::phylo_kmer::score_type _log_threshold;
        const size_t _keep_at_most;
        const double _keep_factor;
        const bool _exists;
    };
}

#endif //RAPPAS_PLACE_PLACE_H
