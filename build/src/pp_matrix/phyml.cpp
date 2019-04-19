#include <iostream>
#include <csv-parser/csv.h>
#include <numeric>
#include "phyml.h"
#include "proba_matrix.h"

using std::string;
using std::cout, std::endl;

class phyml_output_reader
{
public:
    phyml_output_reader(const string& file_name)
        : _file_name(file_name)
    {}

    phyml_output_reader(const phyml_output_reader&) = delete;
    ~phyml_output_reader() = default;

    proba_matrix read()
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
        catch (io::error::integer_overflow& error)
        {
            throw std::runtime_error("PhyML result parsing error: " + string(error.what()));
        }
    }

private:
    proba_matrix read_matrix()
    {
        proba_matrix matrix;

        io::CSVReader<5,
            io::trim_chars<' '>,
            io::no_quote_escape<'\t'>,
            io::throw_on_overflow,
            io::single_and_empty_line_comment<'.'>> _in(_file_name);
        _in.read_header(io::ignore_extra_column, "NodeLabel", "A", "C", "G", "T");

        auto node_label = proba_matrix::NOT_A_LABEL;
        core::phylo_kmer::score_type a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            if (node_label == proba_matrix::NOT_A_LABEL)
            {
                throw std::runtime_error("Node label value is too big: " + std::to_string(node_label));
            }

            /// log-transform the probabilities
            auto new_row = row { { { a, 0 }, { c, 1 }, { g, 2 }, { t, 3 } } };
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
                matrix[node_label] = proba_matrix::mapped_type{node_label, {std::move(new_row)} };
            }
        }
        return matrix;
    }

private:
    string _file_name;
};

proba_matrix load_phyml_probas(const std::string& file_name)
{
    phyml_output_reader reader(file_name);
    return reader.read();
}

extended_mapping load_extended_mapping(const std::string& file_name)
{
    cout << "Loading a node mapping: " + file_name << endl;
    artree_label_mapping mapping;

    io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(file_name);
    in.read_header(io::ignore_extra_column, "original_id", "extended_name");
    std::string extended_name;
    branch_id original_id;
    while (in.read_row(original_id, extended_name))
    {
        mapping[extended_name] = original_id;
    }
    cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
    return mapping;
}

bool is_number(const std::string& s)
{
    auto it = s.begin();
    while (it != s.end() && std::isdigit(*it))
    {
        ++it;
    }
    return !s.empty() && it == s.end();
}

artree_label_mapping load_artree_mapping(const std::string& file_name)
{
    cout << "Loading a node mapping: " + file_name << endl;
    artree_label_mapping mapping;

    io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(file_name);
    in.read_header(io::ignore_extra_column, "extended_label", "ARtree_label");
    std::string extended_label, artree_label;
    while (in.read_row(extended_label, artree_label))
    {
        if (is_number(artree_label))
        {
            mapping[extended_label] = core::phylo_kmer::branch_type(std::stoul(artree_label));
        }
    }
    cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
    return mapping;
}