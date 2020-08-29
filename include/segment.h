#ifndef SEGMENT_H
#define SEGMENT_H

template <typename KeyType, typename PosType, typename Floating = long double>
class Segment
{
public:
    Segment() = default;

    Segment(KeyType start_key, PosType start_pos, KeyType end_key, Floating slope) : start_key(start_key), start_pos(start_pos), end_key(end_key), slope(slope){};

    inline bool operator<(const Segment &s)
    {
        return start_key < s.start_key;
    }

    inline bool operator<(const KeyType &k)
    {
        return start_key < k;
    }

private:
    KeyType start_key;
    PosType start_pos;
    KeyType end_key;
    Floating slope;
};

#endif