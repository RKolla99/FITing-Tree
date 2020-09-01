#ifndef FIT_H
#define FIT_H

#include <cstddef>
#include <cassert>
#include <vector>

#include "segment.h"
#include "piecewise_linear_model.h"
#include "stx/btree.h"

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) : (x) + (error))
#define SUB_ERR(x, error) ((x) <= (error) ? 0 : ((x) - (error)))

template <typename KeyType, uint64_t Error = 64, typename Floating = long double>
class FitingTree
{
    static_assert(Error > 0);

private:
    struct ApproxPos
    {
        uint64_t pos;
        uint64_t hi;
        uint64_t lo;
    };

    size_t n;
    KeyType first_key;
    std::vector<Segment<KeyType, uint64_t>> segments;
    stx::btree<KeyType,
               Segment<KeyType, uint64_t>,
               std::pair<KeyType, Segment<KeyType, uint64_t>>,
               std::greater<KeyType>,
               stx::btree_default_map_traits<KeyType, Segment<KeyType, uint64_t>>,
               false,
               std::allocator<std::pair<KeyType, Segment<KeyType, uint64_t>>>,
               false>
        fiting_tree;

public:
    FitingTree() = default;

    explicit FitingTree(const std::vector<KeyType> &data) : FitingTree(data.begin(), data.end()) {}

    template <typename RandomIt>
    FitingTree(RandomIt first, RandomIt last)
        : n(std::distance(first, last)), first_key(*first), segments(), fiting_tree()
    {
        assert(std::is_sorted(first, last));

        if (n == 0)
            return;

        using pair_type = typename std::pair<KeyType, uint64_t>;
        using tree_pair_type = typename std::pair<KeyType, Segment<KeyType, uint64_t>>;

        std::vector<tree_pair_type> formatted_segments;
        auto error_value = Error;
        size_t num_segments;

        auto in_fun = [this, first](auto i) { return pair_type(first[i], i); };
        auto out_fun = [this](auto segment) { segments.emplace_back(segment); };
        num_segments = get_all_segments(n, error_value, in_fun, out_fun);

        formatted_segments.reserve(num_segments);
        for (auto it = segments.rbegin(); it != segments.rend(); ++it)
        {
            formatted_segments.emplace_back(it->get_start_key(), *it);
        }

        fiting_tree.bulk_load(formatted_segments.begin(), formatted_segments.end());
    }

    ApproxPos get_approx_pos(const KeyType &key)
    {
        if (n == 0)
            return {0, 0, 0};

        auto it = fiting_tree.lower_bound(key);
        if (it == fiting_tree.end())
        {
            return {0, Error, 0};
        }
        else
        {
            KeyType start_key = it->second.get_start_key();
            auto [slope, intercept] = it->second.get_slope_intercept();

            auto pos = (key - start_key) * slope + intercept;

            if (pos - Error > n)
                return {n - 1, n, n - 1};
            return {pos, ADD_ERR(pos, Error, n), SUB_ERR(pos, Error)};
        }
    }
};

#endif