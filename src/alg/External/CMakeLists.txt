add_subdirectory(poisson_disk)

add_library(External)

target_link_libraries(External
    possion_disk)
