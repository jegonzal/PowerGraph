project(GraphLab)

add_graphlab_executable(coloring coloring.cpp)
add_graphlab_executable(factorized_coloring factorized_coloring.cpp)
copy_files(*.tsv)

