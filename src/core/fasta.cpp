#include "fasta.h"
#include <stdexcept>
#include <iostream>
#include <iterator>

using std::vector;
using std::string;
using std::cout, std::endl;
using std::move;

fasta::fasta(string&& header, string&& sequence)
    : _header(move(header))
    , _sequence(move(sequence))
{}

string fasta::get_header() const
{
    return _header;
}

string fasta::get_sequence() const
{
    return _sequence;
}

string fasta::get_sequence_without_gaps() const
{
    throw std::runtime_error("fasta::get_sequence_without_gaps is not supported yet.");
}

//------------------------------------------------------------------------------------
namespace bio = boost::iostreams;

vector<fasta> read_fasta(const string& file_name)
{
    cout << "Loading alignment: " + file_name << endl;
    bio::mapped_file_source mmap(file_name);
    bio::stream<bio::mapped_file_source> is(mmap, std::ios::in);

    vector<fasta> fasta_records;
    string line, header, sequence;
    while (std::getline(is, line))
    {
        if (boost::starts_with(line, ">"))
        {
            if (sequence.size())
            {
                fasta_records.emplace_back(move(header), move(sequence));
                sequence = "";
                sequence.reserve(1024);
            }

            header = move(line);
        }
        else if (!line.empty())
        {
            sequence.append(line);
        }
    }
    fasta_records.emplace_back(move(header), move(sequence));
    return fasta_records;
}