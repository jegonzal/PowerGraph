project(GraphLab)

add_graphlab_executable(shortest_path_cgi shortest_path.cpp)
copy_file(shortest_path.rb)