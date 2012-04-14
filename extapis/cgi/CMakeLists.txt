project(GraphLab)
add_graphlab_executable(dispatcher process.cpp dispatcher.cpp)
copy_files(*.tsv)