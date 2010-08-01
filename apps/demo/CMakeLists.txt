project(GraphLab)

add_executable(demo demo.cpp)
TARGET_LINK_LIBRARIES(demo tcmalloc)
