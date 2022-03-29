#include <fstream>
#include <utility>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include "jplace.h"
#include "place.h"

namespace rappas::io
{
    /// \brief Writes a collection of placements to a .jplace formatted file.
    class jplace_writer
    {
    public:
        explicit jplace_writer(std::string filename) noexcept;
        jplace_writer(const jplace_writer&) = delete;
        jplace_writer(jplace_writer&&) = delete;
        jplace_writer& operator=(const jplace_writer&) = delete;
        jplace_writer& operator=(jplace_writer&&) = delete;
        ~jplace_writer() noexcept = default;

        void write(const std::string& invocation, std::string_view newick_tree,
                   const impl::placed_collection& placed);

    private:
        void _write_metadata(const std::string& invocation);
        void _write_tree(std::string_view newick_tree);
        void _write_version();
        void _write_fields();
        void _write_placements(const impl::placed_collection& placed);
        void _write_placement(const impl::placed_sequence& placed_seq);

        template<class Collection>
        void _write_named_multiplicity(const Collection& seq_headers);

        std::string _filename;
        rapidjson::StringBuffer _buffer;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> _writer;
    };
}

using namespace rappas::io;

jplace_writer::jplace_writer(std::string filename) noexcept
    : _filename{std::move( filename )}, _buffer{ }, _writer{ _buffer }
{}

void jplace_writer::_write_metadata(const std::string& invocation)
{
    _writer.Key("metadata");
    _writer.StartObject();
    _writer.Key("invocation");
    _writer.String(invocation.c_str());
    _writer.EndObject();
}

void jplace_writer::_write_tree(std::string_view newick_tree)
{
    _writer.Key("tree");
    _writer.String(newick_tree.data());
}

void jplace_writer::_write_version()
{
    _writer.Key("version");
    _writer.Uint(3);
}

void jplace_writer::_write_fields()
{
    _writer.Key("fields");
    _writer.StartArray();
    _writer.String("edge_num");
    _writer.String("likelihood");
    _writer.String("like_weight_ratio");
    _writer.String("distal_length");
    _writer.String("pendant_length");
    _writer.EndArray();
}

void jplace_writer::_write_placements(const impl::placed_collection& placed)
{
    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
    _writer.Key("placements");
    _writer.StartArray();
    for (const auto& placed_seq : placed.placed_seqs)
    {
        _writer.StartObject();
        _write_placement(placed_seq);

        const auto seq_headers = placed.sequence_map.at(placed_seq.sequence);
        _write_named_multiplicity(seq_headers);
        _writer.EndObject();
    }
    _writer.EndArray();
}

void jplace_writer::_write_placement(const impl::placed_sequence& placed_seq)
{
    _writer.Key("p");
    _writer.StartArray();
    for (const auto& [branch, score, weight_ratio, count, distal_length, pendant_length] : placed_seq.placements)
    {
        _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
        _writer.StartArray();
        _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatSingleLineArray);
        _writer.Uint(branch);
        _writer.Double(score);
        _writer.Double((double)weight_ratio);
        _writer.Double(distal_length);
        _writer.Double(pendant_length);
        _writer.EndArray();
        (void)count;
    }
    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
    _writer.EndArray();
}

template<class Collection>
void jplace_writer::_write_named_multiplicity(const Collection& seq_headers)
{
    _writer.Key("nm");
    _writer.StartArray();
    for (auto header : seq_headers)
    {
        _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
        _writer.StartArray();
        _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatSingleLineArray);
        _writer.String(header.data());
        _writer.Uint(1);
        _writer.EndArray();
    }
    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
    _writer.EndArray();
}

void jplace_writer::write(const std::string& invocation, std::string_view newick_tree,
    const impl::placed_collection& placed)
{
    /// We use SetFormatOptions to control where to put newlines in the output file.
    /// Possible options are: kFormatSingleLineArray, kFormatDefault
    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatSingleLineArray);

    _writer.StartObject();
    _write_metadata(invocation);
    _write_tree(newick_tree);
    _write_version();
    _write_fields();
    _write_placements(placed);
    _writer.EndObject();

    std::ofstream out(_filename);
    out << _buffer.GetString();
}

void rappas::io::write_jplace(const std::string& filename, const std::string& invocation,
                              std::string_view newick_tree, const impl::placed_collection& placed)
{
    jplace_writer writer(filename);
    writer.write(invocation, newick_tree, placed);
}