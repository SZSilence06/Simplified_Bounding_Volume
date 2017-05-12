project(testCuda)

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    list(APPEND CUDA_NVCC_FLAGS -std=c++11 --expt-relaxed-constexpr -arch=sm_35)
    include_directories(${CUDA_INCLUDE_DIRS})
endif()

file(GLOB headers *.h)
file(GLOB sources *.cu)

cuda_add_executable(testCuda ${headers} ${sources})
