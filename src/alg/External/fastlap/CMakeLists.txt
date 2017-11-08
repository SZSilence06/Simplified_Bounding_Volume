file(GLOB headers *.h)
file(GLOB sources *.c)

add_library(fastlap ${headers} ${sources})
