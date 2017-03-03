
include_directories(".")

file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_subdirectory(mesh)

add_library(wkylib ${headers} ${sources})

target_link_libraries(wkylib
    mesh
    )


