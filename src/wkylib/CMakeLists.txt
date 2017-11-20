cmake_minimum_required(VERSION 2.6)

project(wkylib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

include_directories(".")

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_subdirectory(mesh)
add_subdirectory(Geometry)
add_subdirectory(Cuda)
#add_subdirectory(Test)

add_library(wkylib ${headers} ${sources})

