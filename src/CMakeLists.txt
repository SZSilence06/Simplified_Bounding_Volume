cmake_minimum_required(VERSION 2.8.10)

project(sbvgen)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

find_package(Boost COMPONENTS system filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

include_directories(".")

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_subdirectory(alg)
add_subdirectory(wkylib)

add_executable(sbvgen ${headers} ${sources})

target_link_libraries(sbvgen
    alg
    wkylib
    jtf-mesh
    ${Boost_LIBRARIES}
    )


