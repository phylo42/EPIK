#include "ar.h"
#include "../utils/file_io.h"
#include <iostream>
#include <memory>

using std::string;
using std::cout, std::endl;

proba_matrix::proba_matrix()
{}

void proba_matrix::_add_branch_entry(const string& branch_label, const branch_entry_t& branch_entry)
{
}

class phyml_result_parser
{
public:
    explicit phyml_result_parser(const string& file_name);
    phyml_result_parser(const phyml_result_parser&) = delete;
    ~phyml_result_parser() = default;

    void parse(const string& data);

    proba_matrix get_matrix() const;

private:
    proba_matrix _matrix;
};

phyml_result_parser::phyml_result_parser(const string& file_name)
{}

void phyml_result_parser::parse(const string& data)
{

}

proba_matrix phyml_result_parser::get_matrix() const
{
    return _matrix;
}

proba_matrix load_phyml_probas(size_t seq_alphabet_size, const string& file_name)
{
    phyml_result_parser parser(file_name);
    cout << "Loading PhyML results: " + file_name << endl;

    std::unique_ptr<string_reader> reader = make_mapped_string_reader(file_name);
    if (reader->good())
    {
        string line;
        while (reader->get_line(line))
        {
            parser.parse(line);
            //std::cout << line << std::endl;
        }
    }
    else
    {
        throw std::runtime_error("Cannot open file: " + file_name);
    }
    return parser.get_matrix();
}
