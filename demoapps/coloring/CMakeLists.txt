project(GraphLab)

add_graphlab_executable(coloring coloring.cpp)
copy_files(*.tsv)

