project(GraphLab)


add_distributed_graphlab_executable(dist_als dist_als.cpp)
add_distributed_graphlab_executable(dist_fscope_als dist_als.cpp)
set_target_properties(dist_fscope_als PROPERTIES COMPILE_FLAGS -DFSCOPE)

add_distributed_graphlab_executable(dist_sync_als dist_als.cpp)
set_target_properties(dist_sync_als PROPERTIES COMPILE_FLAGS -DSYNC)


copy_files(*.job)
