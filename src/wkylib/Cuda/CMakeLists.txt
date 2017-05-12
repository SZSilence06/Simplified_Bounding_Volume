find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    list(APPEND CUDA_NVCC_FLAGS -std=c++11 -arch=sm_35)
    include_directories(${CUDA_INCLUDE_DIRS})
endif()

file(GLOB headers *.h)
file(GLOB sources *.cu)

cuda_add_library(wkylib-cuda
    ${headers}
    ${sources})
