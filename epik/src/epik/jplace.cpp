#include <fstream>
#include <utility>
#include <stdexcept>
#include <epik/jplace.h>
#include <epik/place.h>

using namespace epik::io;

jplace_writer::jplace_writer(const std::string& filename,
                             const std::string& invocation,
                             std::string_view newick_tree)
    : _filename(filename), _out(filename), _buffer(),
    _writer(_buffer), _invocation(invocation), _tree(newick_tree)
{
    if (_out.bad())
    {
        throw std::runtime_error("Could not create file " + filename);
    }
}

jplace_writer& jplace_writer::operator<<(const impl::placed_collection& placed)
{
    for (const auto& placed_seq : placed.placed_seqs)
    {
        _writer.StartObject();
        _write_placement(placed_seq);

        const auto seq_headers = placed.sequence_map.at(placed_seq.sequence);
        _write_named_multiplicity(seq_headers);
        _writer.EndObject();
    }

    _out.open(_filename, std::ios_base::app);
    _out << _buffer.GetString();
    _buffer.Clear();
    _out.close();
    return *this;
}

void jplace_writer::start()
{
    /// We use SetFormatOptions to control where to put newlines in the output file.
    /// Possible options are: kFormatSingleLineArray, kFormatDefault
    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatSingleLineArray);

    _writer.StartObject();
    _write_metadata(_invocation);
    _write_tree(_tree);
    _write_version();
    _write_fields();

    _writer.SetFormatOptions(rapidjson::PrettyFormatOptions::kFormatDefault);
    _writer.Key("placements");
    _writer.StartArray();

    _out << _buffer.GetString();
    _buffer.Clear();
    _out.close();
}

void jplace_writer::end()
{
    _writer.EndArray();
    _writer.EndObject();

    _out.open(_filename, std::ios_base::app);
    _out << _buffer.GetString();
    _buffer.Clear();
}

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
