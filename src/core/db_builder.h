#pragma  once

#include <string>
#include <memory>
#include "return.h"
#include "../cli/command_line.h"
#include "seq_traits.h"

class alignment;
class seq_traits;

class db_builder
{
public:
    db_builder(const std::string &tree_file,
               const std::string &working_directory,
               size_t kmer_size,
               const seq_traits &traits,
               double reduction_ratio);

    return_code_t run();

private:
    std::string _tree_file;
    std::string _working_directory;

    size_t _kmer_size;
    seq_traits _seq_traits;
    double _reduction_ratio;
};
