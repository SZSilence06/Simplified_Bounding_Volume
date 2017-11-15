file(GLOB headers *.h *.hpp)
file(GLOB sources *.cpp)

add_library(fmm3d-common
    ${headers}
    ${sources}
)
