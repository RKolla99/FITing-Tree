#ifndef FIT_H
#define FIT_H

#include <cstddef>
#include <vector>

#include "segment.h"
#include "piecewise_linear_model.h"
#include "stx/btree.h"

#define ADD_ERR(x, error, size) ((x) + (error) >= (size) ? (size) : (x) + (error))
#define SUB_ERR(x, error) ((x) <= (error) ? 0 : ((x) - (error)))

template <typename KeyType, size_t Error = 64, typename Floating = long double>
class FitingTree
{
    static_assert(Error > 0);

private:
    size_t n;
    KeyType first_key;
    std::vector<Segment<KeyType, size_t>> segments;
};

#endif