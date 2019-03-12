#pragma once

#include <string>
#include <vector>


class proba_matrix
{
    friend proba_matrix load_phyml_probas(size_t seq_alphabet_size);
    proba_matrix();

public:
    using pos_probs_t = std::vector<float>;
    using branch_probs_t = std::vector<pos_probs_t>;
    using branch_entry_t = std::pair<std::string, branch_probs_t>;

public:
    ~proba_matrix() = default;

private:
    void _add_branch_entry(const std::string& branch_label, const branch_entry_t& branch_entry);

private:
    std::vector<branch_entry_t> _data;
};

proba_matrix load_phyml_probas(size_t seq_alphabet_size);
