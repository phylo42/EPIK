#include <iostream>
#include <csv.h>
#include <numeric>
#include <range/v3/all.hpp>
#include "phyml.h"
#include "proba_matrix.h"

using namespace ranges::v3;

using std::string;
using std::cout, std::endl;

class phyml_output_reader
{
public:
    phyml_output_reader(const string& file_name, const seq_traits& traits)
        : _file_name(file_name)
        , _traits(traits)
    {}

    phyml_output_reader(const phyml_output_reader&) = delete;
    ~phyml_output_reader() = default;

    proba_matrix read()
    {
        try
        {
            cout << "Loading PhyML results: " + _file_name << endl;

            auto matrix = read_matrix();
            cout << "Loaded " << matrix.num_branches() << " matrices of size [" << matrix.num_sites()
                 << " x " << matrix.num_variants() << "]." << endl << endl;

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
        score_t a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            if (node_label == proba_matrix::NOT_A_LABEL)
            {
                throw std::runtime_error("Node label value is too big: " + std::to_string(node_label));
            }

            /// log-transform the probabilities
            auto new_row = row { {a, 0}, {c, 1}, {g, 2}, {t, 3} };
            auto log = [](const proba_pair& p) { return proba_pair{std::log10(p.score), p.index}; };
            std::transform(begin(new_row), end(new_row), begin(new_row), log);

            // sort them
            auto compare = [](const proba_pair& p1, const proba_pair& p2) { return p1.score > p2.score; };
            std::sort(begin(new_row), end(new_row), compare);

            /// insert
            auto it = matrix.find(node_label);
            if (it != end(matrix))
            {
                it->second.push_back(std::move(new_row));
            }
            else
            {
                matrix[node_label] = proba_matrix::mapped_type{node_label, {std::move(new_row)}, dna_seq_traits };
            }
        }
        return matrix;
    }

private:
    string _file_name;
    seq_traits _traits;
};

proba_matrix load_phyml_probas(const std::string& file_name, const seq_traits& traits)
{
    phyml_output_reader reader(file_name, traits);
    return reader.read();
}

node_mapping load_node_mapping(const std::string& file_name)
{
    cout << "Loading a node mapping: " + file_name << endl;
    node_mapping mapping;

    io::CSVReader<2, io::trim_chars<' '>, io::no_quote_escape<'\t'>> in(file_name);
    in.read_header(io::ignore_extra_column, "extended_label", "ARtree_label");
    std::string extended_label, artree_label;
    while (in.read_row(extended_label, artree_label))
    {
        mapping.from_phyml[artree_label] = extended_label;
        mapping.to_phyml[extended_label] = artree_label;
    }
    cout << "Loaded " << mapping.from_phyml.size() << " mapped ids." << endl << endl;
    return mapping;
}