project(GraphLab)



add_graphlab_executable(shortest_path shortest_path.cpp)
copy_files(*.tsv)


