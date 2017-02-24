
file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_library(wkylib_mesh ${headers} ${sources})

