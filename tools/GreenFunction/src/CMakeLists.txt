cmake_minimum_required(VERSION 2.6)

project(GreenFunction)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

find_package(OpenCV)

include_directories("../../../src/")
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(WKY_LIBRARY wkylibd)
    link_directories("../../../build-src-qt5-Clone of Debug/lib/")
elseif(${CMAKE_BUILD_TYPE} MATCHES Release)
    set(WKY_LIBRARY wkylib)
    link_directories("../../../build-src-qt5-Clone of Release/lib/")
endif()


file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_executable(GreenFunction ${headers} ${sources})

target_link_libraries(GreenFunction
    ${OpenCV_LIBS}
    ${WKY_LIBRARY})


