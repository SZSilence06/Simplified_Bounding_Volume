file(GLOB headers *.h)
file(GLOB sources *.cpp)

find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})

add_library(alg ${headers} ${sources})

target_link_libraries(alg
    ${CGAL_LIBRARY}
    ${GMP_LIBRARIES}
    )
