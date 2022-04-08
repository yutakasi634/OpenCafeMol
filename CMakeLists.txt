cmake_minimum_required(VERSION 3.13)
project(open_aicg2plus)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          17)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(src)
