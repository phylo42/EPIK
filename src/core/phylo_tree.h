#pragma once

#include <string>
#include <vector>

class newick_parser;

namespace _impl
{
    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend newick_parser;
    public:
        phylo_node();

        explicit phylo_node(int id, const std::string& label, float branch_length,
                            const std::vector<phylo_node*>& children, phylo_node* parent, bool is_fake);

        phylo_node(const phylo_node& other) = default;
        ~phylo_node() noexcept;

        phylo_node& operator=(const phylo_node&) = delete;

        /// WARNING: this operator only checks for the id and label fields
        bool operator==(const phylo_node& rhs) const noexcept;
        bool operator!=(const phylo_node& rhs) const noexcept;

    private:
        /// Clean node and fill with the default values. Used in the default constructor
        void _clean();

        void add_children(phylo_node* node);

    private:
        int _id;
        std::string _label;
        float _branch_length;
        std::vector<phylo_node*> _children;
        phylo_node* _parent;

        /// Is not used yet. TODO: support this field
        bool _is_fake;
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
    phylo_tree& operator==(const phylo_tree&) = delete;
    phylo_tree&& operator==(phylo_tree&&) = delete;

public:
    ~phylo_tree() noexcept;

    size_t get_node_count() const;

private:
    _impl::phylo_node* _root;
    size_t _node_count;
};

/// \brief A factory function to construct a phylo_tree from a newick file
phylo_tree load_newick(const std::string& file_name);