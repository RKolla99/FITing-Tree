#ifndef BUF_FIT_H
#define BUF_FIT_H

#include <cstddef>
#include <cassert>
#include <vector>
#include <map>

#include "buffered_segment.h"
#include "piecewise_linear_model.h"
#include "stx/btree.h"

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) : (x) + (error))
#define SUB_ERR(x, error) ((x) <= (error) ? 0 : ((x) - (error)))

template <typename KeyType, typename PosType, uint64_t Error = 64, uint64_t BufferSize = 32, typename Floating = long double>
class BufferedFitingTree
{
    static_assert(Error > 0);
    static_assert(BufferSize > 0);
    static_assert(Error > BufferSize);

    class BufferedFitingTreeIterator;

private:
    size_t n;
    KeyType start_key;
    std::vector<BufferedSegment<KeyType, PosType>> segments;
    stx::btree<KeyType,
               BufferedSegment<KeyType, PosType>,
               std::pair<KeyType, BufferedSegment<KeyType, PosType>>,
               std::greater<KeyType>,
               stx::btree_default_map_traits<KeyType, BufferedSegment<KeyType, PosType>>,
               false,
               std::allocator<std::pair<KeyType, BufferedSegment<KeyType, PosType>>>,
               false>
        buffered_fiting_tree;

public:
    static constexpr uint64_t error_value = Error;
    static constexpr uint64_t seg_error = Error - BufferSize;
    static constexpr uint64_t buffer_size = BufferSize;

    using iterator = BufferedFitingTreeIterator;
    using pair_type = typename std::pair<KeyType, PosType>;
    using tree_pair_type = typename std::pair<KeyType, BufferedSegment<KeyType, PosType>>;

    BufferedFitingTree() = default;

    explicit BufferedFitingTree(const std::vector<KeyType> &data) : BufferedFitingTree(data.begin(), data.end()) {}

    template <typename RandomIt>
    BufferedFitingTree(RandomIt first, RandomIt last)
        : n(std::distance(first, last)), start_key(*first), segments(), buffered_fiting_tree()
    {
        assert(std::is_sorted(first, last));

        if (n == 0)
            return;

        std::vector<tree_pair_type> formatted_segments;
        size_t num_segments;

        auto in_fun = [this, first](auto i) { return pair_type(first[i], i); };
        auto out_fun = [this](auto segment) { segments.emplace_back(segment); };
        num_segments = get_all_segments_buffered(n, seg_error, buffer_size, in_fun, out_fun);

        formatted_segments.reserve(num_segments);
        for (auto it = segments.rbegin(); it != segments.rend(); ++it)
        {
            formatted_segments.emplace_back(it->get_start_key(), *it);
        }

        buffered_fiting_tree.bulk_load(formatted_segments.begin(), formatted_segments.end());
    }

    iterator find(const KeyType &key) const
    {
        if (n == 0)
            return end();

        auto it = buffered_fiting_tree.lower_bound(key);
        if (it == buffered_fiting_tree.end())
            return end();

        KeyType start_key = it->second.get_start_key();
        auto [slope, intercept] = it->second.get_slope_intercept();
        auto pos = (key - start_key) * slope;
        auto lower_bound = SUB_ERR(pos, error_value);
        auto upper_bound = ADD_ERR(pos, error_value, it.data().size());

        auto lower_bound_it = it.data().begin();
        auto upper_bound_it = it.data().begin();
        for (; lower_bound > 0; --lower_bound)
            ++lower_bound_it;

        for (; upper_bound > 0; --upper_bound)
            ++upper_bound_it;

        auto segment_it = std::lower_bound(lower_bound_it, upper_bound_it, key);
        if (segment_it != it.data().end() && segment_it->key() == key)
        {
            return segment_it->deleted() ? end() : iterator(this, it, segment_it);
        }

        return end();
    }

    iterator lower_bound(const KeyType &key)
    {
        if (n == 0)
            return end();

        auto it = buffered_fiting_tree.lower_bound(key);
        if (it == buffered_fiting_tree.end())
            return begin();

        KeyType start_key = it.data().get_start_key();
        auto [slope, intercept] = it.data().get_slope_intercept();
        auto pos = (key - start_key) * slope;
        auto lower_bound = SUB_ERR(pos, error_value);
        auto upper_bound = ADD_ERR(pos, error_value, it.data().size());

        auto lower_bound_it = it.data().begin();
        auto upper_bound_it = it.data().begin();
        for (; lower_bound > 0; --lower_bound)
            ++lower_bound_it;

        for (; upper_bound > 0; --upper_bound)
            ++upper_bound_it;

        auto segment_it = std::lower_bound(lower_bound_it, upper_bound_it, key);
        if (segment_it == it.data().end())
        {
            --it;
            if (it == buffered_fiting_tree.rend())
                return end();
            segment_it = it.data().begin();
        }

        while (segment_it->deleted())
        {
            ++segment_it;
            if (segment_it == it.data().end())
            {
                --it;
                if (it == buffered_fiting_tree.rend())
                    return end();
                segment_it = it.data().begin();
            }
        }

        return iterator(this, it, segment_it);
    }

    void insert(const KeyType &key, const PosType &pos)
    {
        if (find(key) != end())
            return;

        auto it = buffered_fiting_tree.lower_bound(key);
        if (it == buffered_fiting_tree.end())
            it = buffered_fiting_tree.rbegin();

        if (!it.data().insert_buffer(key, pos))
        {
            std::vector<pair_type> merged_keys;
            merged_keys = it.data().merge_buffer(key, pos);
            auto merged_keys_it = merged_keys.begin();

            std::vector<BufferedSegment<KeyType, PosType>> new_segments;
            std::vector<tree_pair_type> formatted_segments;
            auto in_fun = [this, merged_keys_it](auto i) { return pair_type(merged_keys_it[i].first, merged_keys_it[i].second); };
            auto out_fun = [&new_segments](auto segment) { new_segments.emplace_back(segment); };

            size_t num_segments = get_all_segments_buffered(merged_keys.size(), seg_error, buffer_size, in_fun, out_fun);
            formatted_segments.reserve(num_segments);
            for (auto it = new_segments.begin(); it != new_segments.end(); ++it)
            {
                formatted_segments.emplace_back(it->get_start_key(), *it);
            }

            auto segment_it = formatted_segments.begin();
            auto current_segment = it.data();
            current_segment = segment_it->second;

            ++segment_it;
            buffered_fiting_tree.insert(segment_it, formatted_segments.end());
        }
    }

    void erase(const KeyType &key)
    {
        auto it = find(key);
        if (it == end())
            return;

        it->set_deleted();
    }

    iterator begin() const
    {
        if (n == 0)
            return end();

        return iterator(this, buffered_fiting_tree.rbegin(), buffered_fiting_tree.rbegin().data().begin());
    }

    iterator end() const
    {
        return iterator(this, buffered_fiting_tree.rend(), buffered_fiting_tree.rbegin().data().end());
    }
};

template <typename K, typename P, uint64_t Error, uint64_t BufferSize, typename Floating>
class BufferedFitingTree<K, P, Error, BufferSize, Floating>::BufferedFitingTreeIterator
{
    friend class BufferedFitingTree;

    using pair_type = typename std::pair<K, P>;
    using segment_iterator = typename BufferedSegment<K, P, Floating>::BufferedSegmentIterator;
    using tree_iterator = typename stx::btree<K,
                                              BufferedSegment<K, P>,
                                              std::pair<K, BufferedSegment<K, P>>,
                                              std::greater<K>,
                                              stx::btree_default_map_traits<K, BufferedSegment<K, P>>,
                                              false,
                                              std::allocator<std::pair<K, BufferedSegment<K, P>>>,
                                              false>::const_reverse_iterator;
    using segment_type = BufferedSegment<K, P, Floating>;
    using buffered_fiting_tree_type = BufferedFitingTree<K, P, Error, BufferSize, Floating>;

    const buffered_fiting_tree_type *super;
    segment_iterator segment_it;
    tree_iterator tree_it;

    void advance_iterator()
    {
        if (tree_it == super->buffered_fiting_tree.rend())
        {
            *this = super->end();
            return;
        }

        ++segment_it;
        if (segment_it == tree_it.data().end())
        {
            ++tree_it;
            if (tree_it == super->buffered_fiting_tree.rend())
            {
                *this = super->end();
                return;
            }
            segment_it = tree_it.data().begin();
        }
    }

    BufferedFitingTreeIterator() = default;

    BufferedFitingTreeIterator(const buffered_fiting_tree_type *super, tree_iterator tree_it, segment_iterator segment_it)
        : super(super), tree_it(tree_it), segment_it(segment_it){};

    BufferedFitingTreeIterator(const BufferedFitingTreeIterator &copy)
        : super(copy.super), tree_it(copy.tree_it), segment_it(segment_it){};

public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const typename BufferedSegment<K, P, Floating>::DataItem;
    using difference_type = const size_t;
    using pointer = const typename BufferedSegment<K, P, Floating>::DataItem *;
    using reference = const typename BufferedSegment<K, P, Floating>::DataItem &;

    BufferedFitingTreeIterator &operator++()
    {
        advance_iterator();
        return *this;
    }

    BufferedFitingTreeIterator operator++(int)
    {
        BufferedFitingTreeIterator i(*this);
        advance_iterator();
        return i;
    }

    reference operator*() const { return *segment_it; }
    pointer operator->() const { return &(*segment_it); }
    bool operator==(const BufferedFitingTreeIterator &rhs) { return ((segment_it == rhs.segment_it) && (tree_it == rhs.tree_it)); }
    bool operator!=(const BufferedFitingTreeIterator &rhs) { return ((segment_it != rhs.segment_it) || (tree_it != rhs.tree_it)); }
};

#endif