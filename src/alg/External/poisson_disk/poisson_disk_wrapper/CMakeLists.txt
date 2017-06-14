file(GLOB headers *.h)
file(GLOB sources *.cpp)

include_directories("../vcglib")

add_library(possion_disk_wrapper ${headers} ${sources})