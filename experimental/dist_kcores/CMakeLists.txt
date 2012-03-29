project(GraphLab)
add_distributed_graphlab_executable(dist_kcores dist_kcores.cpp)
copy_files(*.job)
