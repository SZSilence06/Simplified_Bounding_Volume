add_subdirectory(poisson_disk)
add_subdirectory(fastlap)

add_library(External)

target_link_libraries(External
    possion_disk
    fastlap
    c
    l)
