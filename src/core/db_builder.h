#pragma  once

#include <string>
#include <memory>
#include "return.h"
#include "seq_traits.h"

class alignment;
class seq_traits;

class db_builder
{
public:
    db_builder(const std::string& working_directory,
               const std::string& ar_probabilities_file,
               const std::string& tree_file,
               size_t kmer_size,
               const seq_traits& traits);

    return_code_t run();

private:
    std::string _tree_file;
    std::string _working_directory;

    size_t _kmer_size;
    seq_traits _seq_traits;
};
