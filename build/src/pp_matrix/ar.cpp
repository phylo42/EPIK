#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <string>
#include <sstream>
#include <optional>
#include <regex>
#include <boost/algorithm/string/predicate.hpp>
#include <csv-parser/csv.h>
#include <xpas/seq.h>
#include "ar.h"
#include "proba_matrix.h"


using std::string;
using std::cout, std::endl;
using boost::algorithm::contains;

namespace rappas::io
{
    /// Supported tools for ancestral reconstruction
    enum class ar_format
    {
        PHYML,
        RAXML_NG
    };

    /// \brief Reads ancestral reconstruction output
    class ar_reader
    {
    public:
        virtual ~ar_reader() noexcept = default;
        virtual rappas::proba_matrix read() = 0;
    };

    /// \brief Reads a PhyML output into a matrix.
    class phyml_reader : public ar_reader
    {
    public:
        phyml_reader(const string& file_name) noexcept;
        phyml_reader(const phyml_reader&) = delete;
        phyml_reader(phyml_reader&&) = delete;
        phyml_reader& operator=(const phyml_reader&) = delete;
        phyml_reader& operator=(phyml_reader&&) = delete;
        ~phyml_reader() noexcept override = default;

        proba_matrix read() override;

    private:
        proba_matrix read_matrix();

        string _file_name;
    };

    /// \brief Reads RAXML-NG output into a matrix.
    class raxmlng_reader : public ar_reader
    {
    public:
        raxmlng_reader(const string& file_name) noexcept;
        raxmlng_reader(const phyml_reader&) = delete;
        raxmlng_reader(raxmlng_reader&&) = delete;
        raxmlng_reader& operator=(const raxmlng_reader&) = delete;
        raxmlng_reader& operator=(raxmlng_reader&&) = delete;
        ~raxmlng_reader() noexcept override = default;

        proba_matrix read();

    private:
        proba_matrix read_matrix();

        string _file_name;
    };

    phyml_reader::phyml_reader(const string& file_name) noexcept
        : _file_name{ file_name }
    {}

    proba_matrix phyml_reader::read()
    {
        try
        {
            cout << "Loading PhyML results: " + _file_name << "..." << endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            auto matrix = read_matrix();

            cout << "Loaded " << matrix.num_branches() << " matrices of " <<
                 matrix.num_sites() << " rows." << endl;
            cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << endl << endl;
            return matrix;
        }
        catch (::io::error::integer_overflow& error)
        {
            throw std::runtime_error("PhyML result parsing error: " + string(error.what()));
        }
    }

    proba_matrix phyml_reader::read_matrix()
    {
#ifdef SEQ_TYPE_DNA
        std::ios::sync_with_stdio(false);

        proba_matrix matrix;
        std::ifstream infile(_file_name);

        bool is_header = true;

        string line;
        while (std::getline(infile, line))
        {
            if (is_header)
            {
                // if the line starts with 'Site', the header is finished
                if (line.size() > 4 && line.rfind("Site", 0) == 0)
                {
                    is_header = false;
                }
            }
            else
            {
                size_t site = 0;
                std::string node_label;
                xpas::phylo_kmer::score_type a, c, g, t;

                std::istringstream iss(line);
                if (!(iss >> site >> node_label >> a >> c >> g >> t))
                {
                    throw std::runtime_error("Parsing error: could not parse the line " + line);
                }

                /// log-transform the probabilities
                auto new_row = row_type { { { a, 0 }, { c, 1 }, { g, 2 }, { t, 3 } } };
                auto log = [](const proba_pair& p) { return proba_pair{ std::log10(p.score), p.index }; };
                std::transform(begin(new_row), end(new_row), begin(new_row), log);

                // sort them
                auto compare = [](const proba_pair& p1, const proba_pair& p2) { return p1.score > p2.score; };
                std::sort(begin(new_row), end(new_row), compare);

                /// insert
                auto it = matrix.find(node_label);
                if (it != std::end(matrix))
                {
                    it->second.push_back(std::move(new_row));
                }
                else
                {
                    /// WARNING: it is cheap to copy new_row here in the case of DNA.
                    /// It is probably makes sense to move it in the case of amino acids though.
                    matrix[node_label] = proba_matrix::mapped_type{ node_label, { new_row } };
                }
            }
        }
        return matrix;
#elif SEQ_TYPE_AA
        throw std::runtime_error("PhyML for proteins is not supported yet.");
#else
        static_assert(false, """Make sure the sequence type is defined. Supported types:\n"""
                             """SEQ_TYPE_DNA"""
                             """SEQ_TYPE_AA""");
#endif

    }

    raxmlng_reader::raxmlng_reader(const string& file_name) noexcept
        : _file_name{ file_name }
    {}

    proba_matrix raxmlng_reader::read()
    {
        try
        {
            cout << "Loading RAXML-NG results: " + _file_name << "..." << endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            auto matrix = read_matrix();

            cout << "Loaded " << matrix.num_branches() << " matrices of " <<
                 matrix.num_sites() << " rows." << endl;
            cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << endl << endl;
            return matrix;
        }
        catch (::io::error::integer_overflow& error)
        {
            throw std::runtime_error("PhyML result parsing error: " + string(error.what()));
        }
    }

    proba_matrix raxmlng_reader::read_matrix()
    {
        proba_matrix matrix;

#ifdef SEQ_TYPE_DNA
        ::io::CSVReader<5,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_C", "p_G", "p_T");

        std::string node_label;
        xpas::phylo_kmer::score_type a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            /// log-transform the probabilities
            auto new_row = row_type { { { a, 0 }, { c, 1 }, { g, 2 }, { t, 3 } } };

#elif SEQ_TYPE_AA
        ::io::CSVReader<21,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_R", "p_N", "p_D", "p_C", "p_Q", "p_E", "p_G",
                        "p_H", "p_I", "p_L", "p_K", "p_M", "p_F", "p_P", "p_S", "p_T", "p_W", "p_Y", "p_V");

        std::string node_label;
        xpas::phylo_kmer::score_type a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v;
        while (_in.read_row(node_label, a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v))
        {
            /// log-transform the probabilities
            auto new_row = row_type {
                { { a, *xpas::seq_traits::key_to_code('a') },
                    { r, *xpas::seq_traits::key_to_code('r') },
                    { n, *xpas::seq_traits::key_to_code('n')},
                    { d, *xpas::seq_traits::key_to_code('d') },
                    { c, *xpas::seq_traits::key_to_code('c') },
                    { q, *xpas::seq_traits::key_to_code('q') },
                    { e, *xpas::seq_traits::key_to_code('e') },
                    { g, *xpas::seq_traits::key_to_code('g') },
                    { h, *xpas::seq_traits::key_to_code('h') },
                    { i, *xpas::seq_traits::key_to_code('i') },
                    { l, *xpas::seq_traits::key_to_code('l') },
                    { k, *xpas::seq_traits::key_to_code('k') },
                    { m, *xpas::seq_traits::key_to_code('m') },
                    { f, *xpas::seq_traits::key_to_code('f') },
                    { p, *xpas::seq_traits::key_to_code('p') },
                    { s, *xpas::seq_traits::key_to_code('s') },
                    { t, *xpas::seq_traits::key_to_code('t') },
                    { w, *xpas::seq_traits::key_to_code('w') },
                    { y, *xpas::seq_traits::key_to_code('y') },
                    { v, *xpas::seq_traits::key_to_code('v') } }
            };
#else
            static_assert(false, """Make sure the sequence type is defined. Supported types:\n"""
                             """SEQ_TYPE_DNA"""
                             """SEQ_TYPE_AA""");
#endif
            auto log = [](const proba_pair& p) { return proba_pair{ std::log10(p.score), p.index }; };
            std::transform(begin(new_row), end(new_row), begin(new_row), log);

            // sort them
            auto compare = [](const proba_pair& p1, const proba_pair& p2) { return p1.score > p2.score; };
            std::sort(begin(new_row), end(new_row), compare);

            /// insert into the matrix
            auto it = matrix.find(node_label);
            if (it != std::end(matrix))
            {
                it->second.push_back(std::move(new_row));
            }
            else
            {
                /// WARNING: it is cheap to copy new_row here in the case of DNA.
                /// It probably makes sense to move it in the case of amino acids though.
                matrix[node_label] = proba_matrix::mapped_type{ node_label, { new_row } };
            }
        }
        return matrix;

    }

    ar_format parse_format(const std::string& file_name)
    {
        if (contains(file_name, "phyml"))
        {
            return ar_format::PHYML;
        }
        else if (contains(file_name, "raxml"))
        {
            return ar_format::RAXML_NG;
        }
        throw std::runtime_error("Unsupported ancestral reconstruction result format: " + file_name);
    }

    std::unique_ptr<ar_reader> make_reader(ar_format format, const string& filename)
    {
        if (format == ar_format::PHYML)
        {
            return std::make_unique<phyml_reader>(filename);
        }
        else if (format == ar_format::RAXML_NG)
        {
            return std::make_unique<raxmlng_reader>(filename);
        }
        else
        {
            throw std::runtime_error("Unsupported ancestral reconstruction output format.");
        }
    }

    proba_matrix load_ar(const string& file_name)
    {
        const auto format = parse_format(file_name);
        auto reader = make_reader(format, file_name);
        return reader->read();
    }

    /// \brief Reads a "extended_tree_node_mapping.tsv" file produced by the old RAPPAS.
    extended_mapping load_extended_mapping(const string& file_name)
    {
        cout << "Loading a node mapping: " + file_name << endl;
        extended_mapping mapping;

        ::io::CSVReader<2, ::io::trim_chars<' '>, ::io::no_quote_escape<'\t'>> in(file_name);
        in.read_header(::io::ignore_extra_column, "original_id", "extended_name");
        std::string extended_name;
        branch_type original_id = xpas::phylo_kmer::na_branch;
        while (in.read_row(original_id, extended_name))
        {
            mapping[extended_name] = original_id;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }

    std::optional<xpas::phylo_kmer::branch_type> extract_number(const std::string& s)
    {
        /// take the first sequence of digits from the string
        string output = std::regex_replace(s,
            std::regex("[^0-9]*([0-9]+).*"),
            string("$1")
        );

        /// return it if there was any
        if (output.size() > 0)
        {
            const auto casted = static_cast<xpas::phylo_kmer::branch_type>(std::stoul(output));
            return { casted };
        }
        else
        {
            return std::nullopt;
        }
    }

    /// \brief Returns if the input string is a number
    bool is_number(const std::string& s)
    {
        auto it = s.begin();
        while (it != s.end() && std::isdigit(*it))
        {
            ++it;
        }
        return !s.empty() && it == s.end();
    }

    /// \brief Reads a "ARtree_id_mapping.tsv" file produced by the old RAPPAS.
    artree_label_mapping load_artree_mapping(const std::string& file_name)
    {
        cout << "Loading a node mapping: " + file_name << endl;
        artree_label_mapping mapping;

        ::io::CSVReader<2, ::io::trim_chars<' '>, ::io::no_quote_escape<'\t'>> in(file_name);
        in.read_header(::io::ignore_extra_column, "extended_label", "ARtree_label");
        std::string extended_label, artree_label;
        while (in.read_row(extended_label, artree_label))
        {
            mapping[extended_label] = artree_label;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }
}