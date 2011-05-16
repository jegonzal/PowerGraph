project(GraphLab)

if (ITPP-FOUND)

add_graphlab_executable(gabp gabp.cpp)

#add_graphlab_executable(BlockGABP gabp_block.cpp)

#target_link_libraries(BlockGABP itpp)

endif()
