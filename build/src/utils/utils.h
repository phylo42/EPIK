
#ifndef RAPPAS_BUILD_UTILS_H
#define RAPPAS_BUILD_UTILS_H

#include <vector>
#include <algorithm>
#include <cmath>


namespace rappas {
    namespace utils {
        /// Linear interpolation
        template<typename T>
        static inline T lerp(T v0, T v1, T t) {
            return (1 - t) * v0 + t * v1;
        }

        /// Calculate quantiles of an empirical distribution using linear interpolation
        /// Based on solution from here: https://stackoverflow.com/questions/11964552/finding-quartiles
        template<typename T>
        static inline T quantile(const std::vector<T> &data, T prob) {
            if (data.empty()) {
                return {};
            }

            if (data.size() == 1) {
                return data[0];
            }

            if (!std::is_sorted(data.begin(), data.end())) {
                throw std::runtime_error("Quantile error: input vector must be sorted");
            }

            T poi = lerp<T>(-0.5, data.size() - 0.5, prob);

            const auto left = std::max(static_cast<int64_t>(std::floor(poi)), static_cast<int64_t>(0));
            const auto right = std::min(static_cast<int64_t>(std::ceil(poi)),
                                        static_cast<int64_t>(data.size()) - 1);

            return lerp<T>(data.at(left), data.at(right), poi - left);
        }
    }

}

#endif //RAPPAS_BUILD_UTILS_H
