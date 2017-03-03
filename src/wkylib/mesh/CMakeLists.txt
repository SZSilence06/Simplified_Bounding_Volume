
file(GLOB headers *.h)
file(GLOB sources *.cpp)

#include and link geo_sim_sdk
include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

add_library(mesh ${headers} ${sources})
