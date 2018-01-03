file(GLOB headers *.h)
file(GLOB sources *.cpp)
file(GLOB cuda_sources *.cu)

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
    list(APPEND CUDA_NVCC_FLAGS -arch=sm_50 --default-stream per-thread --expt-relaxed-constexpr -Xcudafe \"--diag_suppress=extra_semicolon\" )
    set(CUDA_SEPARABLE_COMPILATION ON)
endif()

cuda_add_library(alg-cuda ${cuda_sources})

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

#find_package(PCL REQUIRED COMPONENTS common io surface features search filters)
#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#add_definitions(${PCL_DEFINITIONS})

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(External)

add_library(alg ${headers} ${sources})

target_link_libraries(alg
    alg-cuda
    wkylib
    poisson_disk
    fastlap
    lapack
    ${CGAL_LIBRARY}
    ${GMP_LIBRARIES}
#    ${PCL_LIBRARIES}
    ${VTK_LIBRARIES})
