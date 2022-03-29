#include <iostream>
#include <string>
#include <chrono>
#include <boost/filesystem.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>
#include <xpas/phylo_tree.h>
#include <xpas/newick.h>
#include <xpas/fasta.h>
#include <rappas/place.h>
#include <rappas/jplace.h>

namespace fs = boost::filesystem;

/// \brief Creates a string with wich the program was executed
std::string make_invocation(int argc, char** argv)
{
    std::vector<std::string> all_args(argv, argv + argc);
    std::string invocation = "";
    for (const auto& arg : all_args)
    {
        invocation += arg + " ";
    }
    return invocation;
}

fs::path make_output_filename(const std::string& input_file, const std::string& output_dir)
{
    return fs::path(output_dir) / fs::path{ "placements_" + fs::path(input_file).filename().string() + ".jplace" };
}

void print_line()
{
    for (size_t j = 0; j < 60; ++j)
    {
        std::cout << '*';
    }
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);

    if (argc < 5)
    {
        std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE OUTPUT_DIRECTORY EXISTS [0/1] NUM_THREADS QUERY_FILE"  << std::endl;
        return 1;
    }

    const auto keep_at_most = size_t{ 7ul };
    const auto keep_factor = double{ 0.01f };

    const auto db_file = std::string{ argv[1] };
    const auto output_dir = std::string{ argv[2] };
    const auto exists = std::stoul(argv[3]);
    const auto num_threads = std::stoul(argv[4]);

    std::cout << "Loading database..." << std::endl;
    const auto db = xpas::load(db_file);
    std::cout << "Database parameters:" << std::endl
              << "\tSequence type: " << db.sequence_type() << std::endl
              << "\tk: " << db.kmer_size() << std::endl
              << "\tomega: " << db.omega() << std::endl
              << "\tPositions loaded: " << (db.positions_loaded() ? "true" : "false") << std::endl << std::endl;
    std::cout << "Loaded a database of " << db.size() << " phylo-kmers. " << std::endl << std::endl;

    const auto tree = xpas::io::parse_newick(db.tree());
    const auto placer = rappas::placer(db, tree, keep_at_most, keep_factor, exists);
    /// Here we transform the tree to .newick by our own to make sure the output format is always the same
    const auto tree_as_newick = xpas::io::to_newick(tree, true);

    for (int i = 5; i < argc; ++i)
    {
        print_line();
        const auto query_file = std::string{ argv[i] };

        auto sequences = std::vector<xpas::seq_record>();
        for (const auto& seq : xpas::io::read_fasta(query_file))
        {
            sequences.push_back(seq);
        }

        /// Place sequences from a file
        std::cout << "Placing " << query_file << "..." << std::endl;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        const auto placed_seqs = placer.place(sequences, num_threads);
        std::cout << "Placed " << sequences.size() << " sequences.\n" << std::flush;
        std::cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - begin).count() << std::endl << std::endl;

        /// Write results
        const auto jplace_filename = make_output_filename(query_file, output_dir).string();
        const auto invocation = make_invocation(argc, argv);
        std::cout << "Writing to file: " << jplace_filename << "...\n" << std::flush;
        rappas::io::write_jplace(jplace_filename, invocation, tree_as_newick, placed_seqs);
        std::cout << std::endl;
    }

    std::cout << "Done." << '\n' << std::flush;
}