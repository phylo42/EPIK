#include <iostream>
#include <fstream>
#include <csv-parser/csv.h>
#include <numeric>
#include <cmath>
#include <string>
#include <sstream>
#include "phyml.h"
#include "proba_matrix.h"

using std::string;
using std::cout, std::endl;

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


namespace rappas
{
    namespace io
    {
        /// \brief Reads a PhyML output into a matrix.
        class phyml_output_reader
        {
        public:
            phyml_output_reader(const string& file_name) noexcept;
            phyml_output_reader(const phyml_output_reader&) = delete;
            phyml_output_reader(phyml_output_reader&&) = delete;
            phyml_output_reader& operator=(const phyml_output_reader&) = delete;
            phyml_output_reader& operator=(phyml_output_reader&&) = delete;
            ~phyml_output_reader() noexcept = default;

            rappas::proba_matrix read();

        private:
            rappas::proba_matrix read_matrix();

            string _file_name;
        };
    }
}

rappas::io::phyml_output_reader::phyml_output_reader(const string& file_name) noexcept
    : _file_name{ file_name }
{}

rappas::proba_matrix rappas::io::phyml_output_reader::read()
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

rappas::proba_matrix rappas::io::phyml_output_reader::read_matrix()
{
    std::ios::sync_with_stdio(false);

    proba_matrix matrix;
    std::ifstream infile(_file_name);

    bool is_header = true;

    std::string line;
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
            auto node_label = proba_matrix::NOT_A_LABEL;
            xpas::phylo_kmer::score_type a, c, g, t;

            std::istringstream iss(line);
            if (!(iss >> site >> node_label >> a >> c >> g >> t))
            {
                throw std::runtime_error("Parsing error: could not parse the line " + line);
            }

            if (node_label == proba_matrix::NOT_A_LABEL)
            {
                throw std::runtime_error("Node label value is too big: " + std::to_string(node_label));
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
                matrix[node_label] = proba_matrix::mapped_type{node_label, { new_row } };
            }
        }
    }
    return matrix;
}


namespace rappas::io
{
    rappas::proba_matrix load_phyml_probas(const std::string& file_name)
    {
        phyml_output_reader reader(file_name);
        return reader.read();
    }

    /// \brief Reads a "extended_tree_node_mapping.tsv" file produced by the old RAPPAS.
    extended_mapping load_extended_mapping(const std::string& file_name)
    {
        cout << "Loading a node mapping: " + file_name << endl;
        extended_mapping mapping;

        ::io::CSVReader<2, ::io::trim_chars<' '>, ::io::no_quote_escape<'\t'>> in(file_name);
        in.read_header(::io::ignore_extra_column, "original_id", "extended_name");
        std::string extended_name;
        branch_type original_id = xpas::phylo_kmer::nan_branch;
        while (in.read_row(original_id, extended_name))
        {
            mapping[extended_name] = original_id;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
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
            if (is_number(artree_label))
            {
                mapping[extended_label] = xpas::phylo_kmer::branch_type(std::stoul(artree_label));
            }
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }
}