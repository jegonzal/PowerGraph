project(GraphLab)

add_distributed_graphlab_executable(dist_delta_sssp dist_delta_sssp.cpp)




copy_files(*.job)
