add_subdirectory(poisson_disk_wrapper)

add_library(possion_disk)

target_link_libraries(possion_disk
    possion_disk_wrapper)