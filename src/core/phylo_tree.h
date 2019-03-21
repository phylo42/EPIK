#pragma once

#include <string>
#include <vector>
#include <algorithm>

class newick_parser;
class phylo_tree;

namespace _impl
{
    template<bool IsConst>
    class phylo_tree_iterator;

    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend newick_parser;
        template<bool IsConst> friend class phylo_tree_iterator;
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
        void set_label(const std::string& label);

        std::vector<phylo_node*> get_children() const;

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

    template <bool flag, class IsTrue, class IsFalse>
    struct choose;

    template <class IsTrue, class IsFalse>
    struct choose<true, IsTrue, IsFalse> {
        typedef IsTrue type;
    };

    template <class IsTrue, class IsFalse>
    struct choose<false, IsTrue, IsFalse> {
        typedef IsFalse type;
    };

    /// \brief A forward access (non-)const iterator for phylo_node objects. Performs a depth-first
    /// search among a subtree of input phylo_node.
    template <bool IsConst>
    class phylo_tree_iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef typename choose<IsConst, const phylo_node&, phylo_node&>::type reference;
        typedef typename choose<IsConst, const phylo_node*, phylo_node*>::type pointer;

    public:
        phylo_tree_iterator()
            : phylo_tree_iterator(nullptr)
        {}

        phylo_tree_iterator(phylo_node* node)
            : _current(node)
        {}

        phylo_tree_iterator(const phylo_tree_iterator& other) = default;
        ~phylo_tree_iterator()
        {}

        phylo_tree_iterator& operator=(const phylo_tree_iterator& rhs)
        {
            if (*this != rhs)
            {
                _current = rhs._current;
            }
            return *this;
        }

        bool operator==(const phylo_tree_iterator& rhs) const
        {
            return _current == rhs._current;
        }

        bool operator!=(const phylo_tree_iterator& rhs) const
        {
            return !(*this == rhs);
        }

        phylo_tree_iterator& operator++()
        {
            /// Go upside down if necessary. We need to know the index of current node in the parent->children
            phylo_node* temp = _current->_parent;
            int idx = _id_in_parent(_current);
            while (idx == -1 && temp)
            {
                temp = _current->_parent;
                idx = _id_in_parent(_current);
            }

            /// the end of the tree
            if (temp == nullptr)
            {
                _current = nullptr;
            }
            /// visit the next sibling
            else if ((size_t)idx + 1 < temp->_children.size())
            {
                _current = temp->_children[idx + 1];
            }
            /// visit the parent
            else
            {
                _current = temp;
            }
            return *this;
        }

        reference operator*()
        {
            return *_current;
        }

        pointer operator->()
        {
            return _current;
        }

    private:
        int _id_in_parent(phylo_node* node) const
        {
            if (node->_parent != nullptr)
            {
                const auto& children = node->_parent->_children;
                const auto it = std::find(begin(children), end(children), node);
                if (it != end(children))
                {
                    return distance(begin(children), it);
                }
            }
            return -1;
        }
    private:
        phylo_node* _current;
    };
}


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
    phylo_tree& operator=(const phylo_tree&) = delete;
    phylo_tree&& operator=(phylo_tree&&) = delete;

public:
    typedef _impl::phylo_tree_iterator<true> const_iterator;
    typedef _impl::phylo_tree_iterator<false> iterator;

public:
    ~phylo_tree() noexcept;

    size_t get_node_count() const;

    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;

private:
    _impl::phylo_node* _root;
    size_t _node_count;
};

/// A factory function to construct a phylo_tree from a newick file
phylo_tree load_newick(const std::string& file_name);

/// \brief Returns a boolean which indicates if a phylogenetic tree node is fake or not.
/// \details This function just parses a node label. By convention, a node is fake if its label
/// ends with "_X0" or "_X1"
bool is_fake(const _impl::phylo_node& node);