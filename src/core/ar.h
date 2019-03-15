#pragma once

#include <string>
#include <vector>

class phyml_result_parser;

/// \brief A posterior probabilities matrix class
class proba_matrix
{
    friend phyml_result_parser;
public:
    using pos_probs_t = std::vector<float>;
    using branch_probs_t = std::vector<pos_probs_t>;

    struct branch_entry_t
    {
        branch_entry_t();
        branch_entry_t(const branch_entry_t&) = default;
        branch_entry_t(branch_entry_t&&) = default;

        void clear();

    public:
        int node_label;
        branch_probs_t values;
    };

    proba_matrix();
    proba_matrix(const proba_matrix&) = delete;
    proba_matrix(proba_matrix&& other);
    ~proba_matrix() = default;

    proba_matrix& operator=(const proba_matrix&) = delete;

    size_t num_branches() const;
    size_t num_sites() const;
    size_t num_variants() const;

private:
    void _add_branch_entry(const branch_entry_t& branch_entry);

private:
    std::vector<branch_entry_t> _data;
};

proba_matrix load_phyml_probas(const std::string& file_name);
