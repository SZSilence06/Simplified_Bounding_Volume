cmake_minimum_required(VERSION 2.6)

project(sbvgen)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -fopenmp -O3")

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


