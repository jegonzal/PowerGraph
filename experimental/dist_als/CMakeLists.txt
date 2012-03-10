project(GraphLab)

if(HAS_EIGEN)
  add_distributed_graphlab_executable(dist_als dist_als.cpp)
  add_distributed_graphlab_executable(dist_fscope_als dist_als.cpp)
  set_target_properties(dist_fscope_als PROPERTIES COMPILE_FLAGS -DFSCOPE)
endif()

copy_files(*.job)
