cmake_minimum_required(VERSION 2.6)

project(CurveConverter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

include_directories("../../../src/")
link_directories("../../../build-src-qt5-Debug/lib/")

add_executable(CurveConverter main.cpp)

target_link_libraries(CurveConverter
   wkylibd)
