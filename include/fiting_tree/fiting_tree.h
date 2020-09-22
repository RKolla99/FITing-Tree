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

/**
 * A novel index structure that leverages porperties about the underlying data distribution to reduce
 * the size of an index. 
 * 
 * The index is built on a sorted sequence of keys. A query returns a struct ApproxPos containing
 * the approximate position of the key and the bounds of the range of size 2*Error where the key to be
 * found exists, if present. In the case of repeated keys, the index finds the position of the first 
 * occurrence.
 * 
 * The @p Error template parameter should be set according to the desired space-time trade-off. A 
 * smaller error value makes the estimation more precise and the range smaller but at the cost of 
 * increased space usage.
 * 
 * @tparam KeyType - The type of the indexed elements
 * @tparam Error - The maximum error allowed in the segmentation process
 * @tparam Floating - The floating-point type to use for storing slopes
*/
template <typename KeyType, uint64_t Error = 64, typename Floating = long double>
class FitingTree
{
    static_assert(Error > 0);

private:
    /**
     * A struct that stores the result of a query to a @ref FITing-Tree, that is, a range [@ref lo, @ref hi)
     * centered around an approximate position @ref pos of the sought key.
     */
    struct ApproxPos
    {
        uint64_t pos; // The approximate position of the key
        uint64_t hi;  // The lower bound of the range where the key can be found
        uint64_t lo;  // The upper bound of the range where the key can be found
    };

    size_t n;                                         // Total number of keys
    KeyType first_key;                                // The smallest key
    std::vector<Segment<KeyType, uint64_t>> segments; // The segments composing the index
    stx::btree<KeyType,
               Segment<KeyType, uint64_t>,
               std::pair<KeyType, Segment<KeyType, uint64_t>>,
               std::greater<KeyType>,
               stx::btree_default_map_traits<KeyType, Segment<KeyType, uint64_t>>,
               false,
               std::allocator<std::pair<KeyType, Segment<KeyType, uint64_t>>>,
               false>
        fiting_tree; // STX B+ Tree containing all the segments

public:
    /**
     * Constructs an empty index.
     */
    FitingTree() = default;

    /**
     * Constructs the index on the given sorted data.
     * @param data the vector of keys, must be sorted
     */
    explicit FitingTree(const std::vector<KeyType> &data) : FitingTree(data.begin(), data.end()) {}

    /**
     * Constructs the index on the sorted data in the range [first, last).
     * @param first, last the range containing the sorted elements to be indexed
     */
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

    /**
     * Returns the approximate position of a key.
     * @param key the value of the element to search for
     * @return a struct with the approximate position
     */
    ApproxPos get_approx_pos(const KeyType &key) const
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

            uint64_t hi = ADD_ERR(pos, Error, n);
            uint64_t lo = SUB_ERR(pos, Error);
            return {(uint64_t)pos, hi, lo};
        }
    }

    /**
     * Returns the number of segments in the last level of the index.
     * @return the number of segments
     */
    size_t get_segments_count() const
    {
        return segments.size();
    }
};

#endif