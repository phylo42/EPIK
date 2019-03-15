#include "ar.h"
#include "../utils/file_io.h"
#include <iostream>
#include <memory>
#include <boost/algorithm/string.hpp>

using std::string;
using std::cout, std::endl;
using std::vector;
using std::pair;

void swap(proba_matrix::branch_entry_t& b1, proba_matrix::branch_entry_t& b2)
{
    std::swap(b1.node_label, b2.node_label);
    std::swap(b1.values, b2.values);
}

proba_matrix::branch_entry_t::branch_entry_t()
    : node_label(-1)
{}

void proba_matrix::branch_entry_t::clear()
{
    proba_matrix::branch_entry_t empty;
    swap(*this, empty);
}

proba_matrix::proba_matrix()
    : _data()
{}

proba_matrix::proba_matrix(proba_matrix&& other)
    : _data(std::move(other._data))
{}

size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return _data.begin()->values.size();
}

size_t proba_matrix::num_variants() const
{
    return _data.begin()->values.begin()->size();
}

void proba_matrix::_add_branch_entry(const branch_entry_t& branch_entry)
{
    _data.push_back(branch_entry);
}

/// \brief A parser class for the results of PhyML ancestral reconstruction.
class phyml_result_parser
{
public:
    explicit phyml_result_parser();
    phyml_result_parser(const phyml_result_parser&) = delete;
    ~phyml_result_parser() = default;
    phyml_result_parser& operator=(const phyml_result_parser&) = default;

    /// Parse an input buffer. This function can be called more than once,
    /// during the buffered reading from disk.
    /// \param data A string variable containing a current buffer to parse.
    void parse(const string& data);

    /// Get the resulting matrix of posterior probabilities. This function
    /// is supposed to be called at the very end of the parsing process,
    /// after all of the `parse(buffer)` calls.
    /// \sa parse
    proba_matrix&& get_matrix();

private:
    void _process_line(const string& line);
    void _finish_branch();

private:
    /// A temporary storage for the resulting probability matrix
    proba_matrix _matrix;

    string _current_data;
    proba_matrix::branch_entry_t _current_branch;
};

phyml_result_parser::phyml_result_parser()
    : _matrix()
    , _current_data("")
{}

void phyml_result_parser::parse(const string& new_data)
{
    _current_data.append(new_data);

    vector<string> lines;
    boost::split(lines, _current_data, boost::is_any_of("\n"));

    size_t last_line_idx = lines.size() - 1;
    if (_current_data.size() && _current_data[_current_data.size() - 1] != '\n')
    {
        _current_data = lines[last_line_idx];
        --last_line_idx;
    }

    for (size_t i = 0; i <= last_line_idx; ++i)
    {
        _process_line(lines[i]);
    }
}

vector<string> _split_line(const string& line)
{
    auto is_space = [](char c) -> bool { return c == ' ' || c == '\t'; };
    auto is_empty = [](const string& s) -> bool { return s.empty(); };
    vector<string> tokens;
    boost::split(tokens, line, is_space);

    // erase-remove idiom to remove all the spacebars and \t from the array of values
    tokens.erase(std::remove_if(tokens.begin(), tokens.end(), is_empty), tokens.end());
    return tokens;
}

pair<int, proba_matrix::pos_probs_t> _parse_line(const vector<string>& tokens)
{
    auto it = tokens.begin();

    // skip the "Site" column
    ++it;

    // parse the "Node Label" column
    int node_label = std::stoi(*it);
    ++it;

    // parse probability values
    proba_matrix::pos_probs_t probs;
    std::transform(it, tokens.end(), std::back_inserter(probs),
                   [](const string& s) -> float { return std::stof(s); });

    return std::make_pair(node_label, probs);
}

void phyml_result_parser::_process_line(const string& line)
{
    // skip empty lines and the lines that are not fully read yet
    // (they will be parsed during the next call of the parse() method
    if (line.empty() || line == "")
    {
        return;
    }

    vector<string> tokens = _split_line(line);

    // if this line starts from a digit in first column
    string& first_token = tokens[0];
    if (!first_token.empty() && std::isdigit(first_token[0]))
    {
        pair<int, proba_matrix::pos_probs_t> values = _parse_line(tokens);

        // check if it is the first line for a new branch
        if (values.first != _current_branch.node_label)
        {
            _finish_branch();
            _current_branch.node_label = values.first;
        }
        _current_branch.values.push_back(values.second);
    }
}

void phyml_result_parser::_finish_branch()
{
    // put the previous branch to a matrix if exists
    if (!_current_branch.values.empty())
    {
        _matrix._add_branch_entry(_current_branch);
    }
    _current_branch.clear();
}

proba_matrix&& phyml_result_parser::get_matrix()
{
    // We can not know by the file format if the input is over or not. So we assume that we need to
    // process the last branch entry before returning the result;
    // 1) This operation is safe if there is nothing to push.
    // 2) phyml_result_parser::get_matrix is not re-entrant anyway, because it moves the answer
    _finish_branch();

    return std::move(_matrix);
}

proba_matrix load_phyml_probas(const string& file_name)
{
    cout << "Loading PhyML results: " + file_name << endl;

    phyml_result_parser parser;
    buffered_reader reader(file_name);
    if (reader.good())
    {
        while (!reader.empty())
        {
            string chunk = reader.read_next_chunk();
            parser.parse(chunk);
        }
    }
    else
    {
        throw std::runtime_error("Cannot open file: " + file_name);
    }
    proba_matrix matrix = parser.get_matrix();
    cout << "Loaded " << matrix.num_branches() << " matrices of size [" << matrix.num_sites()
         << " x " << matrix.num_variants() << "]." << endl << endl;
    return matrix;
}
