project(GraphLab)



add_graphlab_executable(shortest_path shortest_path.cpp)
add_graphlab_executable(delta_shortest_path delta_shortest_path.cpp)
copy_files(*.tsv)


