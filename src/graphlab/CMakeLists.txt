project(GraphLab)

subdirs(
 distributed
 engine
 extern
 graph
 logger
 monitoring
 parallel
 schedulers
 scope
 serialization
 shared_data
 tasks
 util)
  

add_library(graphlab STATIC
  logger/logger.cpp
  factors/table_factor.cpp
  factors/binary_factor.cpp
  factors/unary_factor.cpp
  parallel/pthread_tools.cpp
  serialization/iarchive.cpp
  serialization/oarchive.cpp
  util/timer.cpp
  util/generics/any.cpp
  util/command_line_options.cpp)



target_link_libraries(graphlab graphlab_metis)


INSTALL(TARGETS graphlab ARCHIVE DESTINATION lib)


# if(MPI_FOUND)
#   add_library(graphlab_dist STATIC
#     distributed/distributed_control.cpp
#     distributed/dc_send.cpp
#     distributed/dc_receive.cpp
#     distributed/dc_internal.cpp
#     distributed/distributed_hash_table.cpp
#     distributed/distributed_terminator.cpp
#     distributed/metrics/distributed_metrics.cpp
#   )

#   # TARGET_LINK_LIBRARIES(graphlab_dist ${MPI_LIBRARY}
#   #  ${MPI_EXTRA_LIBRARY} rt)

#   TARGET_LINK_LIBRARIES(graphlab_dist ${MPI_LIBRARY}
#     ${MPI_EXTRA_LIBRARY} )
#   set_property(TARGET graphlab_dist APPEND PROPERTY COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
# endif(MPI_FOUND)
