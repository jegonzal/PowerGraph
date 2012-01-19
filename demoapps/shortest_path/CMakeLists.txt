project(GraphLab)



add_graphlab_executable(shortest_path shortest_path.cpp)
add_graphlab_executable(coloring coloring.cpp)
copy_directory(graph_data)

