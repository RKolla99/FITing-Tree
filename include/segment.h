#ifndef SEGMENT_H
#define SEGMENT_H

template <typename KeyType, typename PosType, typename Floating = long double>
class Segment
{
private:
    KeyType start_key;
    PosType start_pos;
    KeyType end_key;
    Floating slope;

public:
    Segment() = default;

    Segment(KeyType start_key, PosType start_pos, KeyType end_key, Floating slope) : start_key(start_key), start_pos(start_pos), end_key(end_key), slope(slope){};

    KeyType get_start_key()
    {
        return start_key;
    }

    std::pair<long double, long double> get_slope_intercept()
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