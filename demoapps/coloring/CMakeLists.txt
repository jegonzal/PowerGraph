project(GraphLab)

add_graphlab_executable(coloring coloring.cpp)
add_graphlab_executable(coloring_advanced coloring_advanced.cpp)
copy_files(*.tsv)

