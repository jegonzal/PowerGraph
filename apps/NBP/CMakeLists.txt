project(GraphLab)


add_executable(NBP NBP.cpp)
target_link_libraries(NBP tcmalloc)
target_link_libraries(NBP itpp)
