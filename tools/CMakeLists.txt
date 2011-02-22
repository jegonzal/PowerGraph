project(graphlab_tools)

link_libraries(${Boost_LIBRARIES})
link_libraries(${GraphLab_LIBRARIES})

include_directories( . )

add_dist2_executable(graph_coloring graph_coloring.cpp)
add_dist2_executable(graph_part graph_part.cpp)

