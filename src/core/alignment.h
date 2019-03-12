#pragma once

#include <vector>
#include <string>
#include <memory>
#include "seq_traits.h"

class fasta;

class alignment
{
public:
    typedef std::vector<fasta>::iterator iterator;
    typedef std::vector<fasta>::const_iterator const_iterator;
public:
    explicit alignment(std::vector<fasta>&& sequences);

    // capacity
    size_t width() const;
    size_t height() const;

    // iterators
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
private:
    std::vector<fasta> _sequences;
};

//------------------------------------------------------------------------------------
/*
 * load the alignment from a fasta file
 */
alignment load_alignment(const std::string& file_name);

/*
 * save alignment to a fasta file
 */
void save_alignment(const alignment& align, const std::string& file_name);

/*
 * reduces alignment by deleting all columns containing a proportion of gaps
 * (dash or dot) >= the given ratio
 * note that this operation copies the whole alignment
 */
alignment reduce_alignment(const seq_traits &traits, const alignment &align, double reduction_ratio);

//------------------------------------------------------------------------------------
class alignment_validation
{
public:
    explicit alignment_validation(const seq_traits& seq_traits);
    ~alignment_validation() = default;

    void validate(const alignment& align);
private:
    void _check_length(const alignment& align) const;
    void _check_sequence_states(const alignment& align) const;
    void _check_sequence_states(const fasta& sequence) const;
private:
    seq_traits _seq_traits;
};