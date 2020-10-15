cmake_minimum_required(VERSION 3.0.0)
project(triangles VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(triangles main.cpp triangle.hpp triangle.cpp triangle_space.hpp triangle_space.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
