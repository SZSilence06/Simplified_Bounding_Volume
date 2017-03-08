file(GLOB headers *.h)
file(GLOB sources *.cpp)

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

add_library(alg ${headers} ${sources})
