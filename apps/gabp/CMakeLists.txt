project(GraphLab)


add_executable(GABP gabp.cpp)
target_link_libraries(GABP tcmalloc)

#add_executable(BlockGABP gabp_block.cpp)

#target_link_libraries(BlockGABP itpp)
