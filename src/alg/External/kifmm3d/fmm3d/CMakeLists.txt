file(GLOB headers *.h *.hpp)
file(GLOB sources
    kernel3d.cpp knlmat3d.cpp dense3d.cpp let3d.cpp fmm3d.cpp matmgnt3d.cpp fmm3d_setup.cpp fmm3d_eval.cpp fmm3d_check.cpp)

add_library(fmm3d
    ${headers}
    ${sources}
)
