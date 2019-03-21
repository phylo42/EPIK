#include "phyml.h"
#include "../utils/file_io.h"
#include "phylo_tree.h"
#include "proba_matrix.h"

#include <iostream>
#include <memory>
#include <absl/strings/string_view.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>
#include <csv.h>
#include <numeric>

using std::string;
using std::cout, std::endl;
using std::begin, std::end;

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
            cout << "Loading PhyML results: " + _file_name << endl;

            proba_matrix matrix = read_matrix();
            matrix.sort();

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

        proba_matrix::key_t node_label = proba_matrix::NOT_A_LABEL;
        float a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            if (node_label == proba_matrix::NOT_A_LABEL)
            {
                throw std::runtime_error("Node label value is too big: " + std::to_string(node_label));
            }

            auto new_row = std::make_pair<proba_matrix::row_probs_t, proba_matrix::row_pos_t>({a, c, g, t},
                                                                                              {0, 1, 2, 3});
            auto it = matrix.find(node_label);
            if (it != end(matrix))
            {
                it->second.push_back(std::move(new_row));
            }
            else
            {
                matrix[node_label] = {std::move(new_row)};
            }
        }
        return matrix;
    }

private:
    string _file_name;

};

proba_matrix load_phyml_probas(const string& file_name)
{
    phyml_output_reader reader(file_name);
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