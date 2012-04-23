project(GraphLab)

add_distributed_graphlab_executable(dist_coloring dist_coloring.cpp)
# add_distributed_graphlab_executable(dist_sync_coloring dist_coloring.cpp)
add_distributed_graphlab_executable(dist_fscope_coloring dist_coloring.cpp)
set_target_properties(dist_fscope_coloring PROPERTIES COMPILE_FLAGS -DFSCOPE)
#set_target_properties(dist_sync_coloring PROPERTIES COMPILE_FLAGS -DSYNCHRONOUS_ENGINE)

#add_distributed_graphlab_executable(dist_delta_coloring dist_delta_coloring.cpp)




copy_files(*.job)
