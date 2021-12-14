# Introduction

FITing-Tree is a novel index structure that leverages properties about the underlying data distribution to reduce the size of an index. FITing-Tree was proposed by _Galakatos, Alex, et al. "Fiting-tree: A data-aware index structure." Proceedings of the 2019 International Conference on Management of Data. 2019._ You can learn more about the novel index in their [paper](https://dl.acm.org/doi/abs/10.1145/3299869.3319860).

# Using the Index

FITing-Tree is a header-only library. You can clone the repo and start using it.

```bash
git clone https://github.com/RKolla99/FITing-Tree.git
cd FITing-Tree
```

You can copy the `include/fiting_tree` directory to your project's include path

A sample program showing how to index a vector of random integers using the FITing-Tree.

```cpp
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "fiting_tree/fiting_tree.h"

int main() {
    // Generate some random data
    std::vector<int> data(1000000);
    std::generate(data.begin(), data.end(), std::rand);
    data.push_back(42);
    std::sort(data.begin(), data.end());

    // Construct the index
    const int error_value = 32;
    FitingTree<int, error_value> index(data);

    // Query the index
    auto q = 42;
    auto range = index.get_approx_pos(q);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    std::cout << *std::lower_bound(lo, hi, q);

    return 0;
}
```

# Compiling and running the unit tests

You can build the project and run the tests with

```bash
cmake .
cmake --build .
./test/tests
```

# Design

The design has been made to match with [SOSD](https://github.com/learnedsystems/SOSD). The design has been made while referring to the [PGM Index](https://github.com/gvinciguerra/PGM-index) and contains a lot of similarities in the implementation style.

AS mentioned above, the design has been made to match with SOSD and we've integrated the index with the SOSD Benchmark for performance results. The work can be found [here](https://github.com/reddybhargava/SOSD/tree/dev).

# TODO
1. Add updation and deletion features.
2. Commit hash 20d39375e77a5eff40e4b9c8651e28f742593e93 is the latest stable commit
