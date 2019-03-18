#include "ar.h"
#include "../utils/file_io.h"
#include <iostream>
#include <memory>
#include <absl/strings/string_view.h>
#include <absl/strings/str_cat.h>
#include <absl/strings/str_split.h>

using std::string;
using absl::string_view;
using std::cout, std::endl;
using std::vector;
using std::pair;
using std::begin, std::end;

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
    /// during the buffered reading from disk. Buffer is split by a newline
    /// character, and the lines parsed line by line. If the last line does not
    /// end with a newline character, it is considered as incomplete and will
    /// be ignored.
    /// \param data A string view variable containing a non-owning pointer
    /// to a buffer with new data to parse.
    /// \return The number of parsed characters.
    size_t parse(const string_view& new_data);

    void parse_line(const string_view& line);

    /// Get the resulting matrix of posterior probabilities. This function
    /// is supposed to be called at the very end of the parsing process,
    /// after all of the `parse(buffer)` calls.
    /// \sa parse
    proba_matrix&& get_matrix();

private:
    void _finish_branch();

private:
    /// A temporary storage for the resulting probability matrix
    proba_matrix _matrix;
    proba_matrix::branch_entry_t _current_branch;
};

phyml_result_parser::phyml_result_parser()
    : _matrix()
{}

size_t phyml_result_parser::parse(const string_view& new_data)
{
    vector<string_view> lines = absl::StrSplit(new_data, '\n');

    /// if the last token is not a full line with a newline character,
    /// it will remain not parsed
    bool parse_last_line = (new_data.back() == '\n');
    size_t last_line_idx = lines.size() - 1 - (parse_last_line ? 0 : 1);

    for (size_t i = 0; i <= last_line_idx; ++i)
    {
        parse_line(lines[i]);
    }

    size_t byte_parsed = new_data.size();
    if (last_line_idx != lines.size() - 1)
    {
        byte_parsed -= lines.back().size();
    }
    return byte_parsed;
}

void str_trim(string_view& str_view)
{
    const char* ws = " \t\n";
    str_view.remove_prefix(std::min(str_view.find_first_not_of(ws), str_view.size()));
    str_view.remove_suffix((str_view.size() - 1) - std::min(str_view.find_last_not_of(ws), str_view.size() - 1));
}

pair<int, proba_matrix::pos_probs_t> _parse_line(const vector<string_view>& tokens)
{
    auto it = tokens.begin();

    // skip the "Site" column
    ++it;

    // parse the "Node Label" column
    int node_label = 0;
    if (!absl::SimpleAtoi(*it, &node_label))
    {
        throw std::runtime_error("Wrong PhyML input: " + string(*it));
    }
    ++it;

    // parse probability values
    auto stof = [](const string_view& s) -> float
    {
        float v = 0;
        if (!absl::SimpleAtof(s, &v))
        {
            throw std::runtime_error("Wrong PhyML input: " + string(s));
        }
        return v;
    };
    proba_matrix::pos_probs_t probs;
    probs.reserve(tokens.size());
    std::transform(it, tokens.end(), std::back_inserter(probs), stof);
    return std::make_pair(node_label, probs);
}

void phyml_result_parser::parse_line(const string_view& line)
{
    // skip empty lines and the lines that are not fully read yet
    // (they will be parsed during the next call of the parse() method
    if (line.empty() || line == "")
    {
        return;
    }

    vector<string_view> tokens = absl::StrSplit(line, '\t');
    for (auto& token : tokens)
    {
        str_trim(token);
    }

    // if this line starts from a digit in first column
    string_view first_token = tokens[0];
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
    // 2) phyml_result_parser::get_matrix is not re-entrant anyway, because it moves the answer away
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
        string not_parsed_data;
        while (!reader.empty())
        {
            string_view chunk = reader.read_next_chunk();

            /// continue to parse last line from the previous chunk
            if (not_parsed_data.size() > 0)
            {
                auto it = std::find(begin(chunk), end(chunk), '\n');
                if (it != end(chunk))
                {
                    /// we want to copy '\n' too
                    ++it;

                    string_view end_of_last_line(begin(chunk), std::distance(begin(chunk), it));
                    string last_line;
                    absl::StrAppend(&last_line, not_parsed_data, end_of_last_line);
                    parser.parse_line(last_line);
                    chunk.remove_prefix(end_of_last_line.size());
                }
            }

            size_t parsed_size = parser.parse(chunk);
            not_parsed_data = chunk.substr(parsed_size, std::distance(begin(chunk) + parsed_size, end(chunk)));
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
