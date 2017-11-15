project(testCuda)

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    list(APPEND CUDA_NVCC_FLAGS --expt-relaxed-constexpr -arch=sm_35 -D_MWAITXINTRIN_H_INCLUDED -D__STRICT_ANSI__ -D_FORCE_INLINES)
    include_directories(${CUDA_INCLUDE_DIRS})
endif()

file(GLOB headers *.h)
file(GLOB sources *.cu)

cuda_add_executable(testCuda ${headers} ${sources})
