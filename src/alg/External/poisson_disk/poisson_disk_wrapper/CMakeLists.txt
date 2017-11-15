file(GLOB headers *.h)
file(GLOB sources *.cpp)

include_directories("../vcglib")

add_library(poisson_disk ${headers} ${sources})
