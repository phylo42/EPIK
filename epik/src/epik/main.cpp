#include <iostream>
#include <string>
#include <chrono>
#include <future>
#include <sstream>
#include <iomanip>
#include <cctype>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <cxxopts.hpp>
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
template<typename T>
std::string humanize(T num)
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

/// Parse the humanized RAM size to a number.
/// I know that the name of this function is unfortunate.
size_t dehumanize_ram(const std::string& max_ram)
{
    double value;
    char unit = 0;
    std::stringstream ss(max_ram);

    // Parse the numerical part
    ss >> value;
    if (ss.fail())
    {
        throw std::runtime_error("Can't parse max_ram parameter: wrong numerical part");
    }

    // Check if there is a memory unit
    if (!ss.eof())
    {
        ss >> unit;
        if (ss.fail())
        {
            throw std::runtime_error("Can't parse max_ram parameter: wrong unit");
        }
    }

    switch (std::toupper(unit))
    {
        case 0:
            [[fallthrough]];
        case 'B':
            return static_cast<size_t>(value);
        case 'K':
            return static_cast<size_t>(value * 1024);
        case 'M':
            return static_cast<size_t>(value * 1024 * 1024);
        case 'G':
            return static_cast<size_t>(value * 1024 * 1024 * 1024);
        default:
            throw std::runtime_error("Unknown memory unit.");
    }
}

void check_mu(float mu)
{
    if ((mu < 0.0) || (mu > 1.0))
    {
        throw std::runtime_error("Mu has to a value in [0, 1]");
    }
}


int main(int argc, char** argv)
{
    std::ios::sync_with_stdio(false);

    cxxopts::Options options(argv[0], "Evolutionary Placement with Informative K-mers");
    options.add_options()
        ("d,database", "IPK database", cxxopts::value<std::string>())
        ("q,query", "Input query file (.fasta)", cxxopts::value<std::string>())
        ("j,jobs", "Num threads", cxxopts::value<size_t>()->default_value("1"))
        ("batch-size", "Batch size", cxxopts::value<size_t>()->default_value("2000"))
        ("omega", "Determines the threshold value", cxxopts::value<float>()->default_value("1.5"))
        ("mu", "Proportion of the database to load", cxxopts::value<float>()->default_value("1.0"))
        ("max-ram", "Approximate database size to load, MB", cxxopts::value<std::string>())
        ("o,output-dir", "Output directory", cxxopts::value<std::string>())
        ("keep-at-most", "Number of branches to report", cxxopts::value<size_t>()->default_value("7"))
        ("keep-factor", "Minimum LWR to report", cxxopts::value<double>()->default_value("0.01"))
        ("h,help", "Print usage")
        ;

    if (argc == 1)
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    auto parsed_options = options.parse(argc, argv);
    if (parsed_options.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    try
    {
        const auto db_file = parsed_options["database"].as<std::string>();
        const auto query_file = parsed_options["query"].as<std::string>();
        const auto num_threads = parsed_options["jobs"].as<size_t>();
        const auto batch_size = parsed_options["batch-size"].as<size_t>();
        const auto user_omega = parsed_options["omega"].as<float>();
        const auto user_mu = parsed_options["mu"].as<float>();

        const auto keep_at_most = parsed_options["keep-at-most"].as<size_t>();
        const auto keep_factor = parsed_options["keep-factor"].as<double>();
        const auto output_dir = parsed_options["output-dir"].as<std::string>();

        check_mu(user_mu);

        size_t max_entries = std::numeric_limits<size_t>::max();
        if (parsed_options.count("max-ram"))
        {
            const auto max_ram_string = parsed_options["max-ram"].as<std::string>();
            const auto max_ram = dehumanize_ram(max_ram_string);
            max_entries = max_ram / sizeof(i2l::pkdb_value);
            std::cout << "Max-RAM provided: will be loaded not more than "
                      << humanize(max_entries) << " phylo-k-mers." << std::endl;
        }

#ifndef EPIK_OMP
        if (num_threads != 1)
        {
            std::cerr << "EPIK was complied without OpenMP support and can not be run "
                         "in parallel." << std::endl;
            return -2;
        }
#endif
        std::cout << "Loading database with mu=" << user_mu << " and omega="
                  << user_omega << "..." << std::endl;
        const auto db = i2l::load(db_file, user_mu, user_omega, max_entries);
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
        std::cout << "Loaded " << humanize(db.get_num_entries_loaded())
                  << " of " << humanize(db.get_num_entries_total())
                  << " phylo-k-mers. " << std::endl << std::endl;

        const auto tree = i2l::io::parse_newick(db.tree());
        auto placer = epik::placer(db, tree, keep_at_most, keep_factor, num_threads);
        /// Here we transform the tree to .newick by our own to make sure the output format is always the same
        const auto tree_as_newick = i2l::io::to_newick(tree, true);
        const auto jplace_filename = make_output_filename(query_file, output_dir).string();
        const auto invocation = make_invocation(argc, argv);
        const auto total_fasta_size = fs::file_size(query_file);

        auto jplace = epik::io::jplace_writer(jplace_filename, invocation, tree_as_newick);
        jplace.start();

        print_intruction_set();
        std::cout << "Placing " << query_file << "..." << std::endl;

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
            auto ms_diff = (float)(std::chrono::duration_cast<std::chrono::milliseconds>
                (end_batch - begin_batch).count());
            if (ms_diff == 0)
                ms_diff = 1;
            const auto seq_per_second = 1000.0 * (double)batch_size / ms_diff;
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
        std::cout << "Output: " << jplace_filename << std::endl;

        const auto placement_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - begin).count();
        std::cout << "Placement time: " << humanize_time(placement_time)
            << " (" << placement_time << " ms)" << termcolor::reset << std::endl;
        std::cout << "Done." << '\n' << std::flush;
    }
    catch (const std::runtime_error& error)
    {
        std::cerr << "Error: " << error.what() << std::endl;
        return -1;
    }

    return 0;
}