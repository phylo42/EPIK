#ifndef RAPPAS_CPP_PHYML_H
#define RAPPAS_CPP_PHYML_H

#include <string>
#include <unordered_map>
#include "row.h"

namespace rappas
{
    class proba_matrix;

    using artree_label_mapping = std::unordered_map<std::string, branch_type>;
    using extended_mapping = std::unordered_map<std::string, branch_type>;

    namespace io
    {
        /// Load a posterior probability matrix from a file, generated by PhyML
        rappas::proba_matrix load_phyml_probas(const std::string& file_name);

        /// Load a node mapping from file
        extended_mapping load_extended_mapping(const std::string& file_name);
        artree_label_mapping load_artree_mapping(const std::string& file_name);
    }
}



#endif