project(GraphLab)

add_distributed_graphlab_executable(dist_delta_sssp dist_delta_sssp.cpp)

add_distributed_graphlab_executable(dist_sync_sssp dist_delta_sssp.cpp)
set_target_properties(dist_sync_sssp PROPERTIES COMPILE_FLAGS -DSYNCHRONOUS_ENGINE)

copy_files(*.job)
