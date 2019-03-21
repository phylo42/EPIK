#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <memory>
#include <unordered_map>
#include "return.h"
#include "seq_traits.h"
#include "phylo_kmer.h"

class alignment;
class seq_traits;


using phylo_kmer_db = std::unordered_map<phylo_kmer::kmer_value_t, phylo_kmer>;


class db_builder
{
public:
    db_builder(const std::string& working_directory,
               const std::string& ar_probabilities_file,
               const std::string& tree_file,
               const std::string& mapping_file,
               size_t kmer_size,
               const seq_traits& traits);

    return_code_t run();

private:
    std::string _working_directory;
    std::string _ar_probabilities_file;
    std::string _tree_file;
    std::string _mapping_file;

    size_t _kmer_size;
    seq_traits _seq_traits;
    phylo_kmer_db _phylo_kmer_db;
};

#endif