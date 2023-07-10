#include <iostream>
#include <string>
#include <chrono>
#include <future>
#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <i2l/phylo_kmer_db.h>
#include <i2l/serialization.h>
#include <i2l/phylo_tree.h>
#include <i2l/newick.h>
#include <i2l/fasta.h>
#include <epik/place.h>
#include <epik/jplace.h>

/// \brief Creates a string with wich the program was executed
std::string make_invocation(int argc, char** argv)
{
    std::vector<std::string> all_args(argv, argv + argc);
    std::string invocation;
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
bool is_busy(const std::future<R>& f)
{
    return f.valid() && (f.wait_for(std::chrono::milliseconds(0)) != std::future_status::ready);
}

auto time_diff(std::chrono::steady_clock::time_point begin, std::chrono::steady_clock::time_point end)
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}

void print_intruction_set()
{
#if defined(EPIK_SSE)
    std::cout << "Instruction set: SSE" << std::endl;
#elif defined(EPIK_AVX)
    std::cout << "Instruction set: AVX" << std::endl;
#elif defined(EPIK_AVX2)
    std::cout << "Instruction set: AVX2" << std::endl;
#elif defined(EPIK_AVX512)
    std::cout << "Intruction set: AVX-512" << std::endl;
#else
    std::cout << "Instruction set: scalar" << std::endl;
#endif
}

/// Float-to-humanized string for better output
std::string humanize(double num)
{
    std::ostringstream oss;
    oss.precision(1);

    if (num < 1000.0)
    {
        oss << std::fixed << num;
    }
    else if (num < 1000000.0)
    {
        oss << std::fixed << num / 1000.0 << "K";
    }
    else if (num < 1000000000.0)
    {
        oss << std::fixed << num / 1000000.0 << "M";
    }
    else
    {
        oss << std::fixed << num / 1000000000.0 << "B";
    }

    return oss.str();
}

/// Size_t-to-string that translates milliseconds to humanized time
std::string humanize_time(size_t milliseconds)
{
    // Conversion constants
    const size_t msPerSec = 1000;
    const size_t msPerMin = 60 * msPerSec;
    const size_t msPerHour = 60 * msPerMin;
    const size_t msPerDay = 24 * msPerHour;

    // Calculate time components
    size_t days = milliseconds / msPerDay;
    milliseconds %= msPerDay;
    size_t hours = milliseconds / msPerHour;
    milliseconds %= msPerHour;
    size_t minutes = milliseconds / msPerMin;
    milliseconds %= msPerMin;
    size_t seconds = milliseconds / msPerSec;

    std::ostringstream oss;
    if (days > 0)
    {
        oss << days << " day";
        if (days > 1)
        {
            oss << "s";
        }
        oss << ", ";
    }

    if (hours > 0 || days > 0)
    {
        oss << std::setw(2) << std::setfill('0') << hours << ":";
    }

    oss << std::setw(2) << std::setfill('0') << minutes << ":"
        << std::setw(2) << std::setfill('0') << seconds;

    return oss.str();
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
        const auto batch_size = 1000u;
        const auto user_omega = 1.7f;
        const auto user_mu = 0.5f;
        const auto db = i2l::load(db_file, user_mu, user_omega);
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
            const auto total_fasta_size = fs::file_size(query_file);

            auto jplace = epik::io::jplace_writer(jplace_filename, invocation, tree_as_newick);
            jplace.start();

            std::cout << "Placing " << query_file << "..." << std::endl;
            print_intruction_set();

            using namespace indicators;
            ProgressBar bar{
                option::BarWidth{60},
                option::Start{"["},
                option::Fill{"="},
                option::Lead{">"},
                option::Remainder{" "},
                option::End{"]"},
                option::PrefixText{"Placing "},
                option::ForegroundColor{Color::green},
                option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
                option::MaxProgress{total_fasta_size}
            };

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

            size_t num_seq_placed = 0;
            std::unordered_map<i2l::phylo_kmer::key_type,
                               std::unordered_map<std::string_view, int>> kmer_map;

            double average_speed = 0.0;
            size_t num_iterations = 0;

            /// Batch query reading
            auto reader = i2l::io::batch_fasta(query_file, batch_size);
            while (true)
            {
                // Synchronous reading of the next batch to place
                const auto batch = reader.next_batch();
                if (batch.empty())
                {
                    break;
                }

                // Place in parallel
                const auto begin_batch = std::chrono::steady_clock::now();
                const auto placed_batch = placer.place(batch, num_threads);
                const auto end_batch = std::chrono::steady_clock::now();

                // Compute placement speed, sequences per second
                const auto ms_diff = (float)(std::chrono::duration_cast<std::chrono::milliseconds>
                    (end_batch - begin_batch).count());
                const auto seq_per_second = 1000.0 * batch_size / ms_diff;
                average_speed += seq_per_second;

                // Update progress bar
                bar.set_option(option::PrefixText{humanize(seq_per_second) + " seq/s "});
                bar.set_option(option::PostfixText{std::to_string(num_seq_placed) + " / ?"});
                bar.set_progress(reader.bytes_read());

                // Synchronous output to the .jplace file
                jplace << placed_batch;

                num_seq_placed += batch.size();
                ++num_iterations;
            }
            jplace.end();

            average_speed /= (double)num_iterations;
            bar.set_option(option::PrefixText{"Done. "});
            bar.set_option(option::PostfixText{std::to_string(num_seq_placed)});
            bar.set_progress(reader.bytes_read());

            std::cout << std::endl << termcolor::bold << termcolor::white
                      << "Placed " << num_seq_placed << " sequences.\nAverage speed: "
                      << humanize(average_speed) << " seq/s.\n";

            const auto placement_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count();
            std::cout << "Placement time: " << humanize_time(placement_time)
                << " (" << placement_time << " ms)" << termcolor::reset << std::endl;
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