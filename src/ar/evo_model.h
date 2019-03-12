#pragma once

#include <string>
#include <set>

namespace ar
{
    const std::set<std::string> evo_models = { "JC69", "HKY85", "K80", "F81", "TN93", "GTR" };

    struct evo_model
    {
    public:
        evo_model(const std::string& n, float a, int c);

        std::string name;
        float alpha;
        int categories;
    };

    evo_model get_default_model();
}