cmake_minimum_required(VERSION 2.6)

project(PoissonReconstruct)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

find_package(PCL REQUIRED COMPONENTS common io surface features search filters)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_executable(PoissonReconstruct ${headers} ${sources})

target_link_libraries(PoissonReconstruct
   ${PCL_LIBRARIES})
