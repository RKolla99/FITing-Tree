cmake_minimum_required(VERSION 3.12)
project(FitingTree VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include lib/stx-btree-0.9/include)

enable_testing()
add_subdirectory(test)