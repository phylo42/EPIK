#include "newick.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <stack>
#include <tuple>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/tokenizer.hpp>

using std::vector, std::stack;
using std::string;
using std::move;
using std::cout, std::endl;

using namespace _impl;

phylo_node::phylo_node()
{
    _clean();
}

phylo_node::phylo_node(int id, const std::string& label, float branch_length,
                       const std::vector<phylo_node*>& children, phylo_node* parent, bool is_fake)
    : _id(id)
    , _label(label)
    , _branch_length(branch_length)
    , _children(children)
    , _parent(parent)
    , _is_fake(is_fake)
{}

phylo_node::~phylo_node() noexcept
{
    for (auto child : _children)
    {
        delete child;
    }
}

bool phylo_node::operator==(const phylo_node& rhs) const noexcept
{
    return (_id == rhs._id) && (_label == rhs._label);
}

bool phylo_node::operator!=(const phylo_node& rhs) const noexcept
{
    return !operator==(rhs);
}

void phylo_node::_clean()
{
    _id = -1;
    _label = "";
    _branch_length = 0;
    _children.clear();
    _parent = nullptr;
    _is_fake = false;
}

void phylo_node::add_children(phylo_node* node)
{
    _children.push_back(node);
}

phylo_tree::phylo_tree(_impl::phylo_node* root, size_t node_count)
    : _root(root)
    , _node_count(node_count)
{}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}
size_t phylo_tree::get_node_count() const
{
    return _node_count;
}

/*!
 * \brief A class for parsing .newick-formatted files.
   \details This class parses phylogenetic trees in the newick format. It designed to support
   a buffered reading from disk.
*/
class newick_parser
{
public:
    newick_parser();
    newick_parser(const newick_parser&) = delete;
    newick_parser(newick_parser&&) = delete;
    ~newick_parser() = default;

    /// Parse an input buffer. This function can be called more than once,
    /// during the buffered reading from disk.
    /// \param data A string variable containing the current buffer data to parse.
    void parse(const std::string& data);

    phylo_node* get_root() const;
    size_t get_node_count() const;

private:
    /// Parse next symbol of input data.
    /// \param symbol
    void _parse_symbol(char symbol);

    /// \brief Handles a left parenthesis in input data.
    /// \details A left parenthesis indicates that a new node with children should be created.
    /// We will create it later though, during the _handle_right_parenthesis call
    /// because the phylo_node class has no default constructor for some design reasons.
    /// \sa phylo_node::phylo_node, _handle_right_parenthesis
    void _handle_left_parenthesis();

    /// \details The list of children for "current" parent node is over.
    /// The next symbols are referred to the parent node
    void _handle_right_parenthesis();

    /// \details Node delimiter, we create a node from the text content we collected so far
    void _handle_comma();

    /// \details End of file, we take last node as root
    void _handle_semicolon();

    /// \details Keep reading the current node description
    void _handle_text(char symbol);

    void _start_node();
    phylo_node* _finish_node();
    void _parse_node_text();

private:
    stack<phylo_node*> _node_stack;
    phylo_node* _root;
    int _node_index;
    string _node_text;

    bool _parsing_node;
    bool _end_of_file;
};

newick_parser::newick_parser()
    : _root(nullptr)
    , _node_index(-1)
    , _parsing_node(false)
    , _end_of_file(false)
{}

void newick_parser::parse(const std::string& data)
{
    for (char c : data)
    {
        _parse_symbol(c);

        if (_end_of_file)
        {
            break;
        }
    }
}

phylo_node* newick_parser::get_root() const
{
    return _root;
}

size_t newick_parser::get_node_count() const
{
    return (size_t)_node_index + 1;
}

void newick_parser::_parse_symbol(char symbol)
{
    switch (symbol)
    {
        case '(':
            _handle_left_parenthesis();
            break;
        case ')':
            _handle_right_parenthesis();
            break;
        case ',':
            _handle_comma();
            break;
        case ';':
            _handle_semicolon();
            break;
        default:
            _handle_text(symbol);
            break;
    }
}

void newick_parser::_handle_left_parenthesis()
{
    _start_node();
    _parsing_node = false;
}

void newick_parser::_handle_right_parenthesis()
{
    _finish_node();
    _parsing_node = true;
}

void newick_parser::_handle_comma()
{
    _finish_node();
}

void newick_parser::_handle_semicolon()
{
    _root = _finish_node();
    _end_of_file = true;
}

void newick_parser::_handle_text(char symbol)
{
    /// A node can start from a parenthesis (ex. "(A, B)C") or from a label (ex. "A").
    /// The second case is equal to "()A". In this case we need to create a node as soon
    /// as we read the first symbol
    if (!_parsing_node)
    {
        _start_node();
        _parsing_node = true;
    }

    /// Keep reading the node description
    _node_text.push_back(symbol);
}

void newick_parser::_start_node()
{
    ++_node_index;
    phylo_node* parent = _node_stack.empty() ? nullptr : _node_stack.top();
    _node_stack.push(new phylo_node());
    _node_stack.top()->_id = _node_index;
    _node_stack.top()->_parent = parent;
}

phylo_node* newick_parser::_finish_node()
{
    _parse_node_text();

    /// Add the node to its parent's list
    phylo_node* current_node = _node_stack.top();
    _node_stack.pop();
    if (current_node->_parent != nullptr)
    {
        current_node->_parent->add_children(current_node);
    }

    _parsing_node = false;
    return current_node;
}

void newick_parser::_parse_node_text()
{
    phylo_node* current_node = _node_stack.top();

    // the content can be like "node_name:branch_length", ":branch_length", or just ""
    if (!_node_text.empty())
    {
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        tokenizer tokens(_node_text, boost::char_separator<char>(":"));
        auto it = begin(tokens);

        // if no id presented
        if (!boost::starts_with(_node_text, ":"))
        {
            {
                current_node->_label = *(it++);
            }
        }
        current_node->_branch_length = std::stof(*it);
    }

    // the current node is over
    _node_text.clear();
}

namespace bio = boost::iostreams;

phylo_tree load_newick(const string& file_name)
{
    using std::fpos;
    using std::ifstream;

    newick_parser parser;
    cout << "Loading newick: " + file_name << endl;

    bio::mapped_file_source mmap(file_name);
    bio::stream<bio::mapped_file_source> is(mmap, std::ios::in);
    if (is)
    {
        // get length of file
        is.seekg(0, ifstream::end);
        fpos<mbstate_t> file_length = is.tellg();
        is.seekg(0, ifstream::beg);

        const int buffer_size = 4096;
        char buffer[buffer_size];

        fpos<mbstate_t> read = 0;
        while (read < file_length)
        {
            std::streamsize size_to_read = std::min(file_length - read, static_cast<std::streamoff>(buffer_size - 1));
            is.read(buffer, size_to_read);
            buffer[buffer_size - 1] = '\0';
            read += size_to_read;
            parser.parse(buffer);
        }
        is.close();
        return phylo_tree(parser.get_root(), parser.get_node_count());
    }
    throw std::runtime_error("Cannot open file: " + file_name);

}
