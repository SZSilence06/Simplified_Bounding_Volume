
include_directories(".")

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_library(wkylib ${headers} ${sources})


