
file(GLOB subheaders *.h )
file(GLOB subsources *.cpp )

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(headers
    ${headers}
    ${subheaders}
    PARENT_SCOPE
    )

set(sources
    ${sources}
    ${subsources}
    PARENT_SCOPE
    )
