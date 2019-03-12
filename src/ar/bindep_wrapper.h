#pragma once

#include <string>

namespace ar
{
    class evo_model;

    class bindep_wrapper
    {
    public:

        bindep_wrapper(const std::string& ar_binary, const std::string& alignment_file,
                       const std::string& tree_file, const std::string& ar_parameters, const evo_model& model);
        ~bindep_wrapper();

        void run() const;

    private:
        std::string build_phyml_command(const std::string &ar_binary, const std::string &alignment_file,
                const std::string &tree_file, const std::string &ar_parameters,
                const ar::evo_model &model) const;

    private:
        std::string command;
    };
}
