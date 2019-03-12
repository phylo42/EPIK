#include "alignment.h"
#include "fasta.h"

#include <stdexcept>
#include <iostream>
#include <algorithm>

using std::vector;
using std::string;


//------------------------------------------------------------------------------------
// alignment class
alignment::alignment(vector<fasta>&& sequences)
    : _sequences(move(sequences))
{}

size_t alignment::width() const
{
    return _sequences[0].get_sequence().size();
}

size_t alignment::height() const
{
    return _sequences.size();
}

alignment::iterator alignment::begin()
{
    return _sequences.begin();
}

alignment::iterator alignment::end()
{
    return _sequences.end();
}

alignment::const_iterator alignment::begin() const
{
    return _sequences.begin();
}

alignment::const_iterator alignment::end() const
{
    return _sequences.end();
}
//------------------------------------------------------------------------------------
// miscellaneous functions
alignment load_alignment(const string& file_name)
{
    return alignment(read_fasta(file_name));
}

void save_alignment(const alignment& align, const string& file_name)
{
    save_fasta(std::begin(align), std::end(align), file_name);
}

vector<double> calculate_gap_ratio(const seq_traits& traits, const alignment& align)
{
    vector<double> ratios(align.width(), 0.0f);

    for (const fasta& seq_record : align)
    {
        const string& sequence = seq_record.get_sequence();
        for (size_t i = 0; i < sequence.size(); ++i)
        {
            if (traits.is_gap(sequence[i]))
            {
                ratios[i]++;
            }
        }
    }

    for (double& ratio : ratios)
    {
        ratio /= (double)align.height();
    }
    return ratios;
}

alignment reduce_alignment(const seq_traits &traits, const alignment &align, double reduction_ratio)
{
    // figure out which columns should be removed
    vector<double> gap_ratios = calculate_gap_ratio(traits, align);
    vector<bool> pos_to_remove(gap_ratios.size(), false);
    transform(begin(gap_ratios), end(gap_ratios), pos_to_remove.begin(),
              [reduction_ratio](double ratio) -> bool { return ratio >= reduction_ratio; });

    vector<fasta> reduced_sequences;
    for (const auto& seq_record : align)
    {
        string header = seq_record.get_header();
        string sequence = seq_record.get_sequence();

        // erase-remove idiom: remove the positions, marked as 1 in the pos_to_remove vector
        size_t index = 0;
        auto predicate = [&index, pos_to_remove](char) -> bool { return pos_to_remove[index++]; };
        sequence.erase(
            remove_if(sequence.begin(), sequence.end(), predicate),
            sequence.end()
        );

        reduced_sequences.emplace_back(move(header), move(sequence));
    }
    return alignment(move(reduced_sequences));
}

//------------------------------------------------------------------------------------
// alignment_validation class
alignment_validation::alignment_validation(const seq_traits& seq_traits)
    : _seq_traits(seq_traits)
{}

void alignment_validation::validate(const alignment& align)
{
    // test if the sequence lengths are the same
    _check_length(align);

    // test if the sequences have unsupported states
    _check_sequence_states(align);
}

void alignment_validation::_check_length(const alignment& align) const
{
    for (const auto &sequence : align)
    {
        if (sequence.get_sequence().size() != align.width())
        {
            throw std::runtime_error(
                "Error: Sequences in the input alignment do not have same number of sites."
                "This sequence has a length different from 1st sequence: " + sequence.get_header() +
                " (1st is " + std::begin(align)->get_header() + ")");
        }
    }
}

void alignment_validation::_check_sequence_states(const fasta& sequence) const
{
    for (const auto& state : sequence.get_sequence())
    {
        if (_seq_traits.is_ambiguous(state))
        {
            if (!_seq_traits.is_gap(state))
            {
                // TODO: logging
                //std::cout << "Ambiguous state (char='" << state << "') will be considered as a gap during AR." << std::endl;
            }
        }
        else if (!_seq_traits.is_valid(state))
        {
            throw std::runtime_error("Reference alignment contains a non supported state: " + sequence.get_sequence());
        }
    }
}

void alignment_validation::_check_sequence_states(const alignment& align) const
{
    for (const auto& sequence : align)
    {
        _check_sequence_states(sequence);
    }
}
//------------------------------------------------------------------------------------