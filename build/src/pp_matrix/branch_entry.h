#ifndef RAPPAS_CPP_BRANCH_ENTRY_H
#define RAPPAS_CPP_BRANCH_ENTRY_H

#include <vector>
#include "branch_entry_view.h"

class view_iterator;

/// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
class branch_entry final
{
public:
    using const_iterator = view_iterator;
    using vector_type = std::vector<row>;

    explicit branch_entry() noexcept = default;
    branch_entry(branch_id _id, vector_type&& rows);
    branch_entry(const branch_entry&) = delete;
    branch_entry(branch_entry&&) = default;
    branch_entry& operator=(const branch_entry&) = delete;
    branch_entry& operator=(branch_entry&&) = default;
    ~branch_entry() noexcept = default;

    const_iterator begin(uint32_t kmer_size) const;
    const_iterator end() const;

    void push_back(row&& r);

    size_t get_alignment_size() const;
    branch_id get_branch_label() const;

    const proba_pair& at(size_t position, size_t variant) const;

private:
    branch_id _branch_label;
    vector_type _rows;
};

bool operator==(const branch_entry& lhs, const branch_entry& rhs);

class view_iterator
{
public:
    using iterator_category = std::forward_iterator_tag;
    using const_reference = const branch_entry_view&;
    using const_pointer = const branch_entry_view*;

    view_iterator(branch_entry_view view) noexcept;
    view_iterator(const view_iterator& view) = delete;
    view_iterator(view_iterator&& view) = delete;
    view_iterator& operator=(const view_iterator&) = delete;
    view_iterator& operator=(view_iterator&&) = delete;
    ~view_iterator() = default;

    view_iterator& operator++();

    bool operator==(const view_iterator& rhs) const noexcept;
    bool operator!=(const view_iterator& rhs) const noexcept;
    const_reference operator*() const noexcept;
    const_pointer operator->() const noexcept;
private:
    branch_entry_view _view;
};

#endif //RAPPAS_CPP_BRANCH_ENTRY_H
