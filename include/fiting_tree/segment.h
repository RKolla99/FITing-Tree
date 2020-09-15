#ifndef SEGMENT_H
#define SEGMENT_H

/**
 * The Segment type represents a segment created during segmentation process of the data.
 * The segments are created using the Shrinking Cone Algorithm.
 * 
 * @tparam KeyType - The type of the key to be indexed
 * @tparam PosType - The type of the positions (usually an unsigned integer type)
 * @tparam Floating - The type used to represent floating-point numbers for slope
*/
template <typename KeyType, typename PosType, typename Floating = long double>
class Segment
{
private:
    KeyType start_key; // The smallest key in the segment
    PosType start_pos; // The position of the smallest key
    KeyType end_key;   // The largest key in the segment
    Floating slope;    // The slope of the segment

public:
    Segment() = default;

    /**
     * Constructs a new segment
     * @param start_key - The smallest key in the segment
     * @param start_pos - The position of the smallest key
     * @param end_key - The largest key in the segment
     * @param slope - The slope of the segment
     */
    Segment(KeyType start_key, PosType start_pos, KeyType end_key, Floating slope) : start_key(start_key), start_pos(start_pos), end_key(end_key), slope(slope){};

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

    inline bool operator<(const Segment &s)
    {
        return start_key < s.start_key;
    }

    inline bool operator<(const KeyType &k)
    {
        return start_key < k;
    }
};

#endif