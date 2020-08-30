#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "fiting_tree.h"

#include <type_traits>

TEMPLATE_TEST_CASE("Segmentation algorithm", "", float, double, uint32_t, uint64_t)
{
    const auto error = GENERATE(32, 64, 128);
    std::vector<TestType> data(1000000);
    std::mt19937 engine(42);
    using RandomFunction = std::function<TestType()>;

    if constexpr (std::is_floating_point<TestType>())
    {
        RandomFunction lognormal = std::bind(std::lognormal_distribution<TestType>(0, 0.5), engine);
        RandomFunction exponential = std::bind(std::exponential_distribution<TestType>(1.2), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, lognormal, exponential);
        std::generate(data.begin(), data.end(), rand);
    }
    else
    {
        RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<TestType>(0, 10000), engine);
        RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<TestType>(0, 10000000), engine);
        RandomFunction binomial = std::bind(std::binomial_distribution<TestType>(50000), engine);
        RandomFunction geometric = std::bind(std::geometric_distribution<TestType>(0.8), engine);
        auto rand = GENERATE_COPY(as<RandomFunction>{}, uniform_dense, uniform_sparse, binomial, geometric);
        std::generate(data.begin(), data.end(), rand);
    }

    std::sort(data.begin(), data.end());
    auto segments = get_all_segments(data.begin(), data.end(), error);
    auto it = segments.begin();
    auto [slope, intercept] = it->get_slope_intercept();

    for (auto i = 0; i < data.size(); i++)
    {
        if (i != 0 && data[i] == data[i - 1])
            continue;

        if (std::next(it) != segments.end() && std::next(it)->get_start_key() <= data[i])
        {
            ++it;
            std::tie(slope, intercept) = it->get_slope_intercept();
        }

        auto pos = (data[i] - it->get_start_key()) * slope + intercept;
        auto offset = std::fabs(i - pos);
        REQUIRE(offset <= error + 1);
    }
}