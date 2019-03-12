#include "file_io.h"
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;
using std::string;

buffered_reader::~buffered_reader()
{}

using std::fpos;
using std::ifstream;

/// \brief A buffered reader class for memory mapped files.
class mapped_buffered_reader : public buffered_reader
{
public:
    mapped_buffered_reader(const string& file_name);
    ~mapped_buffered_reader();

    std::string read_next_chunk() override;
    bool empty() const override;
    bool good() const override;
private:
    fpos<mbstate_t> _get_file_legth();
    void _start_reading();
    void _read_next_chunk();

private:
    using mf_source_t = boost::iostreams::mapped_file_source;
    using stream_t = boost::iostreams::stream<mf_source_t>;

    mf_source_t _msource;
    stream_t _stream;

    bool _started;
    fpos<mbstate_t> _file_length;
    fpos<mbstate_t> _read = 0;

    const int _buffer_size = 4096;
    char* _buffer;
};

mapped_buffered_reader::mapped_buffered_reader(const string& file_name)
    : _msource(file_name)
    , _stream(_msource, std::ios::in)
    , _started(false)
    , _file_length(0)
    , _read(0)
    , _buffer(new char[_buffer_size])
{
    _start_reading();
}

mapped_buffered_reader::~mapped_buffered_reader()
{
    delete[] _buffer;
    _stream.close();
}

fpos<mbstate_t> mapped_buffered_reader::_get_file_legth()
{
    _stream.seekg(0, ifstream::end);
    fpos<mbstate_t> file_length = _stream.tellg();
    _stream.seekg(0, ifstream::beg);
    return file_length;
}

string mapped_buffered_reader::read_next_chunk()
{
    if (!_started)
    {
        _start_reading();
    }
    else
    {
        _read_next_chunk();
    }
    return string(_buffer);
}

bool mapped_buffered_reader::empty() const
{
    if (_file_length == 0)
    {
        return 1;
    }

    return _read == _file_length;
}

bool mapped_buffered_reader::good() const
{
    return (bool)_stream;
}

void mapped_buffered_reader::_start_reading()
{
    _file_length = _get_file_legth();
    _started = true;
}

void mapped_buffered_reader::_read_next_chunk()
{
    _buffer[0] = '\0';

    if (_read < _file_length)
    {
        std::streamsize size_to_read = std::min(_file_length - _read, static_cast<std::streamoff>(_buffer_size - 1));
        _stream.read(_buffer, size_to_read);
        _buffer[_buffer_size - 1] = '\0';
        _read += size_to_read;
    }
}

std::unique_ptr<buffered_reader> make_mapped_buffered_reader(const string& file_name)
{
    return std::make_unique<mapped_buffered_reader>(file_name);
}


string_reader::~string_reader()
{}

/// \brief A string reader class for memory mapped files.
class mapped_string_reader : public string_reader
{
public:
    mapped_string_reader(const string& file_name);
    mapped_string_reader(const mapped_string_reader&) = delete;
    ~mapped_string_reader() = default;

    bool get_line(std::string& line) override;
    bool good() const override;

private:
    using mf_source_t = boost::iostreams::mapped_file_source;
    using stream_t = boost::iostreams::stream<mf_source_t>;

    mf_source_t _msource;
    stream_t _stream;
};

mapped_string_reader::mapped_string_reader(const string& file_name)
    : _msource(file_name)
    , _stream(_msource, std::ios::binary)
{}

bool mapped_string_reader::get_line(std::string& line)
{
    return (bool)std::getline(_stream, line);
}

bool mapped_string_reader::good() const
{
    return (bool)_stream;
}
std::unique_ptr<string_reader> make_mapped_string_reader(const std::string& file_name)
{
    return std::make_unique<mapped_string_reader>(file_name);
}