cmake_minimum_required(VERSION 3.20)
project(hw03)

set(CMAKE_CXX_STANDARD 14)

add_executable(hw03 main.cpp)

find_package(CGAL)
include(${CGAL_USE_FILE})