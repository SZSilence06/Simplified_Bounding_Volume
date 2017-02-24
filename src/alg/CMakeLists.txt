file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_library(alg ${headers} ${sources})
