#include <iostream>
#include <string>
#include <chrono>
#include <future>
#include <boost/filesystem.hpp>
#include <i2l/phylo_kmer_db.h>
#include <i2l/serialization.h>
#include <i2l/phylo_tree.h>
#include <i2l/newick.h>
#include <i2l/fasta.h>
#include <epik/place.h>
#include <epik/jplace.h>

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

template<typename R>
bool is_ready_or_dead(const std::future<R>& f)
{
    return (!f.valid()) || (f.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready);
}

int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);

    if (argc < 4)
    {
        std::cout << "Usage:\n\t" << argv[0] << " DATABASE_FILE OUTPUT_DIRECTORY NUM_THREADS QUERY_FILE" << std::endl;
        return 1;
    }

    try
    {
        const auto keep_at_most = size_t{7ul};
        const auto keep_factor = double{0.01f};

        const auto db_file = std::string{argv[1]};
        const auto output_dir = std::string{argv[2]};
        const auto num_threads = std::stoul(argv[3]);

#ifndef EPIK_OMP
        if (num_threads != 1)
        {
            std::cerr << "EPIK was complied without OpenMP support and can not be run "
                         "in parallel." << std::endl;
            return -2;
        }
#endif

        std::cout << "Loading database..." << std::endl;
        const auto db = i2l::load(db_file);
        if (db.version() < i2l::protocol::EARLIEST_INDEX)
        {
            std::cerr << "The serialization protocol version is too old (v" << db.version() << ").\n"
                      << "Can not use databases built by xpas older than v0.3.2" << std::endl;
            return -1;
        }

        std::cout << "Database parameters:" << std::endl
                  << "\tSequence type: " << db.sequence_type() << std::endl
                  << "\tk: " << db.kmer_size() << std::endl
                  << "\tomega: " << db.omega() << std::endl
                  << "\tPositions loaded: " << (db.positions_loaded() ? "true" : "false") << std::endl << std::endl;
        std::cout << "Loaded a database of " << db.size() << " phylo-kmers. " << std::endl << std::endl;

        const auto tree = i2l::io::parse_newick(db.tree());
        auto placer = epik::placer(db, tree, keep_at_most, keep_factor, num_threads);
        /// Here we transform the tree to .newick by our own to make sure the output format is always the same
        const auto tree_as_newick = i2l::io::to_newick(tree, true);


        for (int i = 4; i < argc; ++i)
        {
            print_line();
            const auto query_file = std::string{argv[i]};

            const auto jplace_filename = make_output_filename(query_file, output_dir).string();
            const auto invocation = make_invocation(argc, argv);

            auto jplace = epik::io::jplace_writer(jplace_filename, invocation, tree_as_newick);
            jplace.start();

            std::cout << "Placing " << query_file << "..." << std::endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            size_t num_seq_placed = 0;

            /// Batch query reading
            auto reader = i2l::io::read_fasta(query_file);
            auto it = reader.begin();
            auto read_batch = [&it, &reader]() {
                const size_t batch_size = 10000;
                auto sequences = std::vector<i2l::seq_record>();
                bool end = (it == reader.end());
                size_t j = 0;
                while (!end)
                {
                    sequences.push_back(*it);
                    ++j;
                    ++it;
                    end = (j >= batch_size) || (it == reader.end());
                }
                return sequences;
            };

            auto write_batch = [&jplace] (const auto& batch) {
                jplace << batch;
            };

            /// Query reading future policy depends on how many threads we are allowed to use.
            /// If only one, read synchronously and place in one thread.
            /// Otherwise, read asynchronously in a separate thread and place with N-1 threads
            const auto async_policy = (num_threads > 1) ? std::launch::async : std::launch::deferred;
            auto future_read = std::async(async_policy, read_batch);
            std::future<void> future_write;

            bool end = false;
            while (!end)
            {
                const auto batch = future_read.get();
                if (batch.empty())
                {
                    end = true;
                }
                future_read = std::async(async_policy, read_batch);

                const auto threads_available = num_threads -
                    (is_ready_or_dead(future_read) ? 0 : 1) - (is_ready_or_dead(future_write) ? 0 : 1);
                std::cout << "Threads: " << threads_available << ". " <<
                    "Reading: " << (is_ready_or_dead(future_read) ? ". " : "BUSY ") <<
                    "Writing: " << (is_ready_or_dead(future_write) ? ". " : "BUSY ") << std::endl;
                const auto placed_batch = placer.place(batch, threads_available);
                if (future_write.valid())
                {
                    future_write.wait();
                }
                const auto async_write_policy = std::launch::async;
                future_write = std::async(async_write_policy, write_batch, placed_batch);

                num_seq_placed += batch.size();
            }
            future_write.wait();
            jplace.end();

            std::cout << "Placed " << num_seq_placed << " sequences.\n" << std::flush;
            std::cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << std::endl << std::endl;
        }

        std::cout << "Done." << '\n' << std::flush;
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << "Exception occurred:\n\t" << error.what() << std::endl;
        return -1;
    }

    return 0;
}