#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "fiting_tree.h"
#include "buffered_fiting_tree.h"

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

TEMPLATE_TEST_CASE_SIG("Fiting-Tree Index", "",
                       ((typename T, size_t E), T, E),
                       (uint32_t, 16), (uint32_t, 32), (uint32_t, 64),
                       (uint64_t, 16), (uint64_t, 32), (uint64_t, 64))
{
    std::vector<T> data(2000000);
    std::mt19937 engine(42);

    using RandomFunction = std::function<T()>;
    RandomFunction uniform_dense = std::bind(std::uniform_int_distribution<T>(0, 10000), engine);
    RandomFunction uniform_sparse = std::bind(std::uniform_int_distribution<T>(0, 10000000), engine);
    RandomFunction binomial = std::bind(std::binomial_distribution<T>(50000), engine);
    RandomFunction geometric = std::bind(std::geometric_distribution<T>(0.8), engine);
    auto rand = GENERATE_COPY(as<RandomFunction>{}, uniform_dense, uniform_sparse, binomial, geometric);

    std::generate(data.begin(), data.end(), rand);
    std::sort(data.begin(), data.end());
    FitingTree<T, E> fiting_tree(data);

    for (auto i = 1; i <= 10000; ++i)
    {
        auto q = data[std::rand() % data.size()];
        auto approx_range = fiting_tree.get_approx_pos(q);
        auto lo = data.begin() + approx_range.lo;
        auto hi = data.begin() + approx_range.hi;
        auto k = std::lower_bound(lo, hi, q);
        REQUIRE(*k == q);
    }

    // Test elements outside range
    auto q = data.back() + 42;
    auto approx_range = fiting_tree.get_approx_pos(q);
    auto lo = data.begin() + approx_range.lo;
    auto hi = data.begin() + approx_range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.end());

    q = 0;
    approx_range = fiting_tree.get_approx_pos(q);
    lo = data.begin() + approx_range.lo;
    hi = data.begin() + approx_range.hi;
    REQUIRE(std::lower_bound(lo, hi, q) == data.begin());
}

TEST_CASE("Buffered Fiting-Tree Iterator")
{
    std::srand(42);
    auto gen = [] { return std::rand() % 1000000000; };

    std::vector<uint32_t> bulk(1000000);
    std::generate(bulk.begin(), bulk.end(), gen);
    std::sort(bulk.begin(), bulk.end());

    BufferedFitingTree<uint32_t, uint32_t> fiting_tree(bulk);

    int i = 0;
    for (auto it = fiting_tree.begin(); it != fiting_tree.end(); ++it)
    {
        REQUIRE(it->key() == bulk[i]);
        i += 1;
    }
}

TEMPLATE_TEST_CASE("Buffered FITing-Tree Index", "", uint32_t, uint64_t)
{
    std::srand(42);
    auto gen = [] { return std::rand() % 1000000000; };

    std::vector<uint32_t> bulk(1000000);
    std::generate(bulk.begin(), bulk.end(), gen);
    std::sort(bulk.begin(), bulk.end());

    BufferedFitingTree<uint32_t, TestType> fiting_tree(bulk);

    for (auto i = 1; i <= 1000; ++i)
    {
        auto q = bulk[std::rand() % bulk.size()];
        auto it = fiting_tree.lower_bound(q);
        REQUIRE(it->key() == q);
    }
}