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
    db_builder(const std::string& ar_binary,
               const std::string& alignment_file,
               const std::string& tree_file,
               const std::string& working_directory,
               const std::string& ar_directory,
               const std::string& ar_parameters,
               size_t kmer_size,
               const ar::evo_model& model,
               const seq_traits& traits,
               double reduction_ratio);
    return_code_t run();

private:
    std::string build_phyml_command() const;

private:
    std::string _ar_binary;
    std::string _alignment_file;
    std::string _tree_file;
    std::string _working_directory;
    std::string _ar_directory;
    std::string _ar_parameters;

    size_t _kmer_size;
    ar::evo_model _model;
    seq_traits _seq_traits;
    double _reduction_ratio;
};
