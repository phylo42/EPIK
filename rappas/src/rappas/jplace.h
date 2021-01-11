#ifndef RAPPAS_PLACE_JPLACE_H
#define RAPPAS_PLACE_JPLACE_H

#include <string>

namespace rappas
{
    namespace impl
    {
        class placed_collection;
    }

    namespace io
    {
        /// \brief Write a collection of placed sequences to a .jplace formatted file
        void write_jplace(const std::string& filename, const std::string& invocation,
                          std::string_view newick_tree, const impl::placed_collection& placed);

    }
}

#endif //RAPPAS_PLACE_JPLACE_H
