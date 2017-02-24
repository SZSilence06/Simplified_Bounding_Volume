cmake_minimum_required(VERSION 2.6)

project(sbvgen)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_subdirectory(alg)
add_subdirectory(wkylib)

add_executable(SBV_APP ${headers} ${sources})

target_link_libraries(SBV_APP
    alg)


