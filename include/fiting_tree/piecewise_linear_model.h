#ifndef PLM_H
#define PLM_H

#include <vector>
#include <cstddef>
#include <type_traits>

#include "segment.h"

template <typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template <typename X, typename Y>
class PiecewiseLinearModel
{
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Slope
    {
        SX dx{};
        SY dy{};

        inline bool operator<(const Slope &p) const
        {
            return dy * p.dx < dx * p.dy;
        }

        inline bool operator>(const Slope &p) const
        {
            return dy * p.dx > dx * p.dy;
        }

        inline bool operator==(const Slope &p) const
        {
            return dy * p.dx == dx * p.dy;
        }

        inline bool operator!=(const Slope &p) const
        {
            return dy * p.dx != dx * p.dy;
        }

        explicit operator long double() const
        {
            return dy / (long double)dx;
        }
    };

    struct Point
    {
        X x{};
        SY y{};

        inline Slope operator-(const Point &p) const
        {
            return {SX(x) - p.x, y - p.y};
        }
    };

    const Y error;
    Point first_point;
    Point last_point;
    Slope lower_slope = {1, 0};
    Slope upper_slope = {0, 1};
    size_t points_in_segment = 0;

public:
    explicit PiecewiseLinearModel(Y error) : error(error)
    {
        if (error < 0)
        {
            throw std::invalid_argument("error can't be less than zero");
        }
    }

    bool add_point(const X &x, const Y &y)
    {
        Point current_point{x, y};
        Point p1{x, SY(y) + error};
        Point p2{x, SY(y) - error};

        if (points_in_segment == 0)
        {
            first_point = current_point;
            last_point = current_point;
            lower_slope = {1, 0};
            upper_slope = {0, 1};
            ++points_in_segment;
            return true;
        }

        if (points_in_segment == 1)
        {
            lower_slope = p2 - first_point;
            upper_slope = p1 - first_point;
            ++points_in_segment;
            last_point = current_point;
            return true;
        }

        Slope slope = current_point - first_point;
        bool outside_lower_slope = slope < lower_slope;
        bool outside_upper_slope = slope > upper_slope;

        if (outside_lower_slope || outside_upper_slope)
        {
            points_in_segment = 0;
            return false;
        }

        if (p1 - first_point < upper_slope)
        {
            upper_slope = p1 - first_point;
        }

        if (p2 - first_point > lower_slope)
        {
            lower_slope = p2 - first_point;
        }

        last_point = current_point;
        ++points_in_segment;
        return true;
    }

    Segment<X, Y> get_segment()
    {
        if (points_in_segment == 1)
            return Segment<X, Y>((X)first_point.x, (Y)first_point.y, (X)last_point.x, 1);
        long double u_slope = (long double)upper_slope;
        long double l_slope = (long double)lower_slope;
        long double slope = (u_slope + l_slope) / 2;
        return Segment<X, Y>(X(first_point.x), Y(first_point.y), X(last_point.x), slope);
    }
};

template <typename Fin, typename Fout>
size_t get_all_segments(size_t n, size_t error, Fin in, Fout out)
{
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t num_segments = 0;
    size_t start = 0;
    auto kv = in(0);

    PiecewiseLinearModel<X, Y> plm(error);
    plm.add_point(kv.first, kv.second);

    for (size_t i = 1; i < n; ++i)
    {
        auto next_kv = in(i);
        if (i != start && next_kv.first == kv.first)
            continue;

        kv = next_kv;
        if (!plm.add_point(kv.first, kv.second))
        {
            out(plm.get_segment());
            start = i;
            --i;
            ++num_segments;
        }
    }

    out(plm.get_segment());
    return ++num_segments;
}

template <typename RandomIterator>
auto get_all_segments(RandomIterator first, RandomIterator last, size_t error)
{
    using key_type = typename RandomIterator::value_type;
    using pair_type = typename std::pair<key_type, size_t>;

    size_t n = std::distance(first, last);
    std::vector<Segment<key_type, size_t>> out;

    auto in_fun = [first](auto i) { return pair_type(first[i], i); };
    auto out_fun = [&out](auto segment) { out.push_back(segment); };
    get_all_segments(n, error, in_fun, out_fun);

    return out;
}

#endif