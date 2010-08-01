project(GraphLab)



add_executable(shooting shooting.cpp)
target_link_libraries(shooting tcmalloc)
