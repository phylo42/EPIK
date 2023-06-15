#ifndef EPIK_JPLACE_H
#define EPIK_JPLACE_H

#include <string>

namespace epik
{
    namespace impl
    {
        struct placed_collection;
    }

    namespace io
    {
        /// \brief Write a collection of placed sequences to a .jplace formatted file
        void write_jplace(const std::string& filename, const std::string& invocation,
                          std::string_view newick_tree, const impl::placed_collection& placed);

    }
}

#endif
