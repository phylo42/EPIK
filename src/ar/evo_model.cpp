#include "evo_model.h"

namespace ar
{
    evo_model::evo_model(const std::string& n, float a, int c)
        : name(n)
        , alpha(a)
        , categories(c)
    {}

    evo_model get_default_model()
    {
        return evo_model("JC69", 1.0f, 4);
    }

}