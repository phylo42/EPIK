#ifndef EPIK_JPLACE_H
#define EPIK_JPLACE_H

#include <string>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>

namespace epik
{
    namespace impl
    {
        class placed_sequence;
        class placed_collection;
    }

    namespace io
    {
        /// \brief Writes a collection of placements to a .jplace formatted file.
        class jplace_writer
        {
        public:
            jplace_writer(const std::string& filename, const std::string& invocation, std::string_view newick_tree);
            jplace_writer(const jplace_writer&) = delete;
            jplace_writer(jplace_writer&&) = delete;
            jplace_writer& operator=(const jplace_writer&) = delete;
            jplace_writer& operator=(jplace_writer&&) = delete;
            ~jplace_writer() noexcept = default;

            jplace_writer& operator<<(const impl::placed_collection& placed);

            void start();
            void end();

        private:
            void _write_metadata(const std::string& invocation);
            void _write_tree(std::string_view newick_tree);
            void _write_version();
            void _write_fields();
            void _write_placements(const impl::placed_collection& placed);
            void _write_placement_batch(const impl::placed_collection& placed);
            void _write_placement(const impl::placed_sequence& placed_seq);

            template<class Collection>
            void _write_named_multiplicity(const Collection& seq_headers);

            std::string _filename;
            std::ofstream _out;
            rapidjson::StringBuffer _buffer;
            rapidjson::PrettyWriter<rapidjson::StringBuffer> _writer;

            std::string _invocation;
            std::string_view _tree;
        };

    }
}

#endif
