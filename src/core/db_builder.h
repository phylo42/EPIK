#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <memory>

#include "return.h"
#include "seq.h"
#include "phylo_kmer_db.h"
#include "pp_matrix/phyml.h"

class alignment;
class branch_entry;
class phylo_tree;
class proba_matrix;

class db_builder
{
public:
    db_builder(std::string working_directory, std::string ar_probabilities_file, std::string tree_file,
               std::string extended_mapping_file, std::string artree_mapping_file, size_t kmer_size);

    return_code run();

private:
    void explore_kmers(const phylo_tree& tree, const proba_matrix& probas);
    size_t explore_branch(const branch_entry& probas, branch_node_t common_branch_label);

    std::string _working_directory;
    std::string _ar_probabilities_file;
    std::string _tree_file;
    std::string _extended_mapping_file;
    std::string _artree_mapping_file;

    size_t _kmer_size;
    phylo_kmer_db _phylo_kmer_db;
    extended_mapping _extended_mapping;
    artree_label_mapping _artree_mapping;
};

#endif