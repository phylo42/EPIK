#include "file_io.h"
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>

namespace bio = boost::iostreams;
using std::string;
using std::fpos;
using std::ifstream;
using absl::string_view;

buffered_reader::buffered_reader(const string& file_name)
    : _msource(file_name)
    , _stream(_msource, std::ios::in)
    , _started(false)
    , _file_length(0)
    , _read(0)
    , _buffer(new char[_buffer_size])
{
    _start_reading();
}

buffered_reader::~buffered_reader()
{
    delete[] _buffer;
    _stream.close();
}

fpos<mbstate_t> buffered_reader::_get_file_legth()
{
    _stream.seekg(0, ifstream::end);
    fpos<mbstate_t> file_length = _stream.tellg();
    _stream.seekg(0, ifstream::beg);
    return file_length;
}

string_view buffered_reader::read_next_chunk()
{
    if (!_started)
    {
        _start_reading();
    }
    else
    {
        _read_next_chunk();
    }
    return string_view(_buffer);
}

bool buffered_reader::empty() const
{
    if (_file_length == 0)
    {
        return 1;
    }

    return _read == _file_length;
}

bool buffered_reader::good() const
{
    return (bool)_stream;
}

void buffered_reader::_start_reading()
{
    _file_length = _get_file_legth();
    _started = true;
}

void buffered_reader::_read_next_chunk()
{
    _buffer[0] = '\0';

    if (_read < _file_length)
    {
        std::streamsize size_to_read = std::min(_file_length - _read, static_cast<std::streamoff>(_buffer_size - 1));
        _stream.read(_buffer, size_to_read);
        _buffer[size_to_read] = '\0';
        _read += size_to_read;
    }
}