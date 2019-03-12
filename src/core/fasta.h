#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string/predicate.hpp>


class fasta
{
public:
    fasta() = delete;
    fasta(std::string&& header, std::string&& sequence);
    fasta(fasta&& other) = default;
    fasta(const fasta&) = delete;
    ~fasta() = default;

    std::string get_header() const;
    std::string get_sequence() const;
    std::string get_sequence_without_gaps() const;
private:
    const std::string _header;
    const std::string _sequence;
};

std::vector<fasta> read_fasta(const std::string& file_name);

// save a stream of fasta records to a file
template <class InputIt>
void save_fasta(InputIt first, InputIt last, const std::string& file_name)
{
    std::cout << "Saving alignment to " << file_name << "..." << std::endl;

    // calculate file size
    size_t size = 0;
    for (auto it = first; it != last; ++it)
    {
        size += it->get_header().size() + it->get_sequence().size() + 2;
    }

    // output itself
    namespace bio = boost::iostreams;
    bio::mapped_file_params params;
    params.path = file_name;
    params.new_file_size = size;
    params.flags = bio::mapped_file::mapmode::readwrite;

    bio::stream<bio::mapped_file_sink> out(params);
    for (auto it = first; it != last; ++it)
    {
        out << it->get_header() << std::endl << it->get_sequence() << std::endl;
    }
}