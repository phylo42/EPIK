#pragma once

#include <string>
#include <vector>

class newick_parser;

class phylo_tree;

namespace _impl
{
    class phylo_tree_iterator;

    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend newick_parser;
        friend phylo_tree_iterator;
        friend phylo_tree;
    public:
        phylo_node();

        explicit phylo_node(int id, const std::string& label, float branch_length,
                            const std::vector<phylo_node*>& children, phylo_node* parent);

        phylo_node(const phylo_node& other) = delete;
        ~phylo_node() noexcept;

        phylo_node& operator=(const phylo_node&) = delete;

        /// WARNING: this operator only checks for the id and label fields
        bool operator==(const phylo_node& rhs) const noexcept;
        bool operator!=(const phylo_node& rhs) const noexcept;

        std::string get_label() const;

    private:
        /// Clean node and fill with the default values. Used in the default constructor
        void _clean();

        void _add_children(phylo_node* node);

    private:
        int _id;
        std::string _label;
        float _branch_length;
        std::vector<phylo_node*> _children;
        phylo_node* _parent;
    };

    /// \brief A constant forward access iterator for phylo_node objects. Performs a depth-first
    /// search among children of a phylo_node.
    class phylo_tree_iterator
    {
    public:
        //typedef int difference_type;
        //typedef A::value_type value_type;
        typedef const phylo_node& reference;
        typedef const phylo_node* pointer;
        typedef std::forward_iterator_tag iterator_category;

    public:
        phylo_tree_iterator();
        phylo_tree_iterator(phylo_node* node);
        phylo_tree_iterator(const phylo_tree_iterator& other);
        //phylo_tree_iterator(const iterator&);
        ~phylo_tree_iterator();

        phylo_tree_iterator& operator=(const phylo_tree_iterator& rhs);
        bool operator==(const phylo_tree_iterator&rhs) const;
        bool operator!=(const phylo_tree_iterator&rhs) const;

        phylo_tree_iterator& operator++();

        reference operator*() const;
        pointer operator->() const;

    private:
        int _id_in_parent(phylo_node* node) const;
    private:
        phylo_node* _current;
    };
}

/// Returns a boolean which indicates if a phylogenetic tree node is fake or not.
/// WARNING: this function just parses a node label. A node is fake if its label
/// ends with "_X0" or "_X1"
bool is_fake(const _impl::phylo_node& node);

/// \brief A phylogenetic tree class
/// \defails phylo_tree is only constructable with a factory function. Also non-copyable and non-movable.
/// \sa load_newick
class phylo_tree
{
private:
    friend phylo_tree load_newick(const std::string& file_name);
    phylo_tree(_impl::phylo_node* root, size_t node_count);
    phylo_tree(phylo_tree&&) = delete;
    phylo_tree(const phylo_tree&) = delete;
    phylo_tree& operator==(const phylo_tree&) = delete;
    phylo_tree&& operator==(phylo_tree&&) = delete;

public:
    typedef _impl::phylo_tree_iterator const_iterator;

public:
    ~phylo_tree() noexcept;

    size_t get_node_count() const;

    const_iterator begin() const;
    const_iterator end() const;

private:
    _impl::phylo_node* _root;
    size_t _node_count;
};

/// A factory function to construct a phylo_tree from a newick file
phylo_tree load_newick(const std::string& file_name);

/// \brief Returns if a phylo_node is fake or not.
/// \details By convention, it checks if node label ends with "_X0" or "_X1".
bool is_fake(const _impl::phylo_node& node);