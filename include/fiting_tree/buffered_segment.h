#ifndef BUF_SEGMENT_H
#define BUF_SEGMENT_H

#include <vector>
#include <map>

/**
 * The BufferedSegment type represents a segment created during segmentation process of the data.
 * The segments are created using the Shrinking Cone Algorithm. It differs from the normal Segment
 * by providing a buffer for every segment to accomodate inserts effeciently.
 * 
 * @tparam KeyType - The type of the key to be indexed
 * @tparam PosType - The type of the positions (usually an unsigned integer type)
 * @tparam Floating - The type used to represent floating-point numbers for slope
*/
template <typename KeyType, typename PosType, typename Floating = long double>
class BufferedSegment
{
public:
    class BufferedSegmentIterator;
    class DataItem;

    using iterator = BufferedSegmentIterator;
    using pair_type = std::pair<KeyType, PosType>;

private:
    KeyType start_key;          // The smallest key in the segment
    PosType start_pos;          // The position of the smallest key
    KeyType end_key;            // The largest key in the segment
    Floating slope;             // The slope of the segment
    std::vector<DataItem> keys; // Stores all the key value pairs present in the segment
    std::map<KeyType,
             DataItem,
             std::less<KeyType>>
        buffer;               // A buffer maintained in sorted order for inserts
    uint64_t buffer_size;     // Current Buffer size
    uint64_t max_buffer_size; // Maximum Buffer size allowed

public:
    BufferedSegment() = default;

    /**
     * Constructs a new segment
     * @param start_key - The smallest key in the segment
     * @param start_pos - The position of the smallest key
     * @param end_key - The largest key in the segment
     * @param slope - The slope of the segment
     */
    BufferedSegment(KeyType start_key, PosType start_pos, KeyType end_key, Floating slope, std::vector<pair_type> p_keys, uint64_t max_buffer_size)
        : start_key(start_key), start_pos(start_pos), end_key(end_key), slope(slope), buffer(), buffer_size(0), max_buffer_size(max_buffer_size)
    {
        keys.reserve(p_keys.size());

        for (auto i = 0; i < p_keys.size(); i++)
        {
            keys.emplace_back(p_keys[i].first, p_keys[i].second);
        }
    }

    /**
     * Returns the smallest key in the segment
     * @return the smallest key 
     */
    KeyType get_start_key() const
    {
        return start_key;
    }

    /**
     * Returns the slope and the intercept of the segment
     * @return a std::pair of [slope, intercept]
     */
    std::pair<long double, long double> get_slope_intercept() const
    {
        return {slope, start_pos};
    }

    bool insert_buffer(const KeyType &key, const PosType &pos)
    {
        if (buffer_size >= max_buffer_size)
            return false;

        buffer.insert({key, DataItem(key, pos)});
        buffer_size += 1;
        return true;
    }

    std::vector<pair_type> merge_buffer(const KeyType &new_key, const PosType &new_pos) const
    {
        std::vector<pair_type> merged_keys;
        merged_keys.reserve(keys.size() + buffer_size + 1);
        bool new_key_added = false;

        auto it = begin();
        while (it != end())
        {
            if (it->deleted())
            {
                ++it;
                continue;
            }

            if (!new_key_added && new_key < it->key())
            {
                merged_keys.emplace_back(new_key, new_pos);
                new_key_added = true;
                continue;
            }

            merged_keys.emplace_back(it->key(), it->pos());
            ++it;
        }

        return merged_keys;
    }

    size_t size() const
    {
        return (keys.size() + buffer_size);
    }

    iterator begin() const
    {
        if (keys.empty())
            return end();

        auto key_it = keys.begin();
        for (auto it = keys.begin(); it != keys.end(); ++it)
        {
            if (it->deleted())
                continue;
            key_it = it;
            break;
        }

        auto buffer_it = buffer.begin();
        for (auto it = buffer.begin(); it != buffer.end(); ++it)
        {
            if (it->second.deleted())
                continue;
            buffer_it = it;
            break;
        }

        return iterator(this, key_it, buffer_it);
    }

    iterator end() const
    {
        return iterator(this, keys.end(), buffer.end());
    }

    inline bool operator<(const BufferedSegment &s)
    {
        return start_key < s.start_key;
    }

    inline bool operator<(const KeyType &k)
    {
        return start_key < k;
    }
};

template <typename K, typename P, typename Floating>
class BufferedSegment<K, P, Floating>::BufferedSegmentIterator
{
    friend class BufferedSegment;

    using item_type = typename BufferedSegment<K, P, Floating>::DataItem;
    using keys_iterator = typename std::vector<item_type>::const_iterator;
    using buffer_iterator = typename std::map<K, item_type>::const_iterator;
    using buffered_segment_type = BufferedSegment<K, P, Floating>;

    const buffered_segment_type *super;
    keys_iterator key_it;
    buffer_iterator buffer_it;

    void advance_iterator()
    {
        if (key_it == super->keys.end() && buffer_it == super->buffer.end())
        {
            *this = super->end();
            return;
        }
        else if (key_it == super->keys.end())
        {
            ++buffer_it;
            return;
        }
        else if (buffer_it == super->buffer.end())
        {
            ++key_it;
            return;
        }
        else
        {
            if (key_it->key() > buffer_it->first)
            {
                ++buffer_it;
                return;
            }
            ++key_it;
        }
    }

public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const DataItem;
    using difference_type = size_t;
    using pointer = const DataItem *;
    using reference = const DataItem &;

    BufferedSegmentIterator() = default;

    BufferedSegmentIterator(const buffered_segment_type *super, keys_iterator key_it, buffer_iterator buffer_it)
        : super(super), key_it(key_it), buffer_it(buffer_it){};

    BufferedSegmentIterator(const BufferedSegmentIterator &copy)
        : super(copy.super), key_it(copy.key_it), buffer_it(copy.buffer_it){};

    BufferedSegmentIterator &operator++()
    {
        advance_iterator();
        return *this;
    }

    BufferedSegmentIterator operator++(int)
    {
        BufferedSegmentIterator i(*this);
        advance_iterator();
        return i;
    }

    reference operator*() const
    {
        if (key_it == super->keys.end() && buffer_it == super->buffer.end())
        {
            return *super->end();
        }
        else if (key_it == super->keys.end())
        {
            return buffer_it->second;
        }
        else if (buffer_it == super->buffer.end())
        {
            return *key_it;
        }
        else
        {
            if (key_it->key() > buffer_it->first)
            {
                return buffer_it->second;
            }
            return *key_it;
        }
    }

    pointer operator->() const
    {
        return &(**this);
    }

    bool operator==(const BufferedSegmentIterator &rhs) const { return ((key_it == rhs.key_it) && (buffer_it == rhs.buffer_it)); }
    bool operator!=(const BufferedSegmentIterator &rhs) const { return ((key_it != rhs.key_it) || (buffer_it != rhs.buffer_it)); }
};

template <typename K, typename P, typename Floating>
class BufferedSegment<K, P, Floating>::DataItem
{
private:
    K first;
    P second;
    bool is_deleted;

public:
    DataItem() = default;
    explicit DataItem(const K &key, const P &pos) : first(key), second(pos), is_deleted(false){};

    bool deleted() const { return is_deleted; }
    void set_deleted() { is_deleted = true; }

    const K &key() const { return first; }
    const P &pos() const { return second; }

    bool operator<(const K &key) const { return first < key; }
};

#endif