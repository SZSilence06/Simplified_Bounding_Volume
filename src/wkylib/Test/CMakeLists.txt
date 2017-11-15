project(wkylib_test)

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_subdirectory(Cuda)
