file(GLOB headers *.h)
file(GLOB sources *.cpp)
file(GLOB cuda_sources *.cu)

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    list(APPEND CUDA_NVCC_FLAGS "-std=c++11")
    include_directories(${CUDA_INCLUDE_DIRS})
endif()

cuda_compile(CUDA_OBJECTS ${cuda_sources})
add_library(alg ${headers} ${sources} ${CUDA_OBJECTS})

target_link_libraries(alg
    wkylib
    ${CGAL_LIBRARY}
    ${GMP_LIBRARIES}
    ${CUDA_LIBRARIES})
