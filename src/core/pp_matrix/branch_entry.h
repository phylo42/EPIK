#ifndef RAPPAS_CPP_BRANCH_ENTRY_H
#define RAPPAS_CPP_BRANCH_ENTRY_H

#include <vector>
#include <core/seq_traits.h>
#include <utils/meta.h>
#include "row.h"
#include "branch_entry_view.h"

class view_iterator;

/// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
class branch_entry final
{
public:
    using const_iterator = view_iterator;

    explicit branch_entry(const seq_traits& traits = dna_seq_traits);
    branch_entry(branch_id _id, std::vector<row>&& rows, const seq_traits& traits);
    branch_entry(const branch_entry&) = delete;
    branch_entry(branch_entry&&) = default;
    branch_entry& operator=(const branch_entry&) = delete;
    branch_entry& operator=(branch_entry&&) = default;
    ~branch_entry() noexcept = default;

    const_iterator begin(size_t kmer_size) const;
    const_iterator end() const;

    void push_back(row&& r);

    size_t get_alignment_size() const;
    size_t get_alphabet_size() const;
    branch_id get_branch_label() const;

    const seq_traits& traits() const;

    const proba_pair& at(size_t position, size_t variant) const;

private:
    branch_id _branch_label;
    std::vector<row> _rows;
    seq_traits _traits;
};

bool operator==(const branch_entry& lhs, const branch_entry& rhs);

class view_iterator
{
public:
    typedef std::forward_iterator_tag iterator_category;
    using reference = const branch_entry_view&;
    using pointer = const branch_entry_view*;

    view_iterator(branch_entry_view view) noexcept;
    view_iterator(const view_iterator& view) = delete;
    view_iterator(view_iterator&& view) = delete;
    view_iterator& operator=(const view_iterator&) = delete;
    view_iterator& operator=(view_iterator&&) = delete;
    ~view_iterator() = default;

    view_iterator& operator++();

    bool operator==(const view_iterator& rhs) const noexcept;
    bool operator!=(const view_iterator& rhs) const noexcept;
    reference operator*();
    pointer operator->();
private:
    branch_entry_view _view;
};

#endif //RAPPAS_CPP_BRANCH_ENTRY_H
