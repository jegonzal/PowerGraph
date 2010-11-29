project(GraphLab)

subdirs(
 distributed
 engine
 extern
 graph
 gpu
 monitoring
 parallel
 schedulers
 scope
 shared_data
 tasks
 util
 serialization
 logger)


add_library(graphlab_util STATIC
  util/timer.cpp
  util/generics/any.cpp
  util/command_line_options.cpp
  logger/logger.cpp
  parallel/pthread_tools.cpp
  schedulers/scheduler_list.cpp
)



if(Sctp-FOUND)
  set(sctp_source rpc/dc_sctp_comm.cpp)
endif()
  
add_library(graphlab_rpc STATIC
  rpc/dc_tcp_comm.cpp
  ${sctp_source}
  rpc/circular_char_buffer.cpp
  rpc/dc_stream_send.cpp
  rpc/dc_stream_receive.cpp
  rpc/dc.cpp
  rpc/reply_increment_counter.cpp
  rpc/dc_services.cpp
  rpc/dc_comm_services.cpp)

if(Sctp-FOUND)
  set_property(TARGET graphlab_rpc APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
endif()

target_link_libraries(graphlab_rpc graphlab_util)


add_library(graphlab STATIC
  factors/table_factor.cpp
  factors/binary_factor.cpp
  factors/unary_factor.cpp
  schedulers/scheduler_options.cpp)
  

if(Sctp-FOUND)
  target_link_libraries(graphlab_rpc sctp)
endif()

target_link_libraries(graphlab graphlab_metis graphlab_rpc graphlab_util)


INSTALL(TARGETS graphlab ARCHIVE DESTINATION lib)


if(MPI_FOUND)
   add_library(graphlab_dist STATIC
     distributed/distributed_control.cpp
     distributed/dc_send.cpp
     distributed/dc_receive.cpp
     distributed/dc_internal.cpp
     distributed/distributed_hash_table.cpp
     distributed/distributed_terminator.cpp
     distributed/metrics/distributed_metrics.cpp
)

 TARGET_LINK_LIBRARIES(graphlab_dist ${MPI_LIBRARY}
  ${MPI_EXTRA_LIBRARY} rt)

   TARGET_LINK_LIBRARIES(graphlab_dist ${MPI_LIBRARY}
     ${MPI_EXTRA_LIBRARY} )
   set_property(TARGET graphlab_dist APPEND PROPERTY COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
 endif(MPI_FOUND)


# PIC versions of the graphlab libraries
get_target_property(graphlab_util_sources graphlab_util SOURCES)
get_property(graphlab_util_flags TARGET graphlab_util PROPERTY COMPILE_FLAGS)
add_library(graphlab_util_pic STATIC ${graphlab_util_sources})
set_target_properties(graphlab_util_pic PROPERTIES 
                      COMPILE_FLAGS "${graphlab_util_flags} -fPIC")



get_target_property(graphlab_sources graphlab SOURCES)
get_property(graphlab_flags TARGET graphlab PROPERTY COMPILE_FLAGS)
add_library(graphlab_pic STATIC ${graphlab_sources})
set_target_properties(graphlab_pic PROPERTIES 
                      COMPILE_FLAGS "${graphlab_flags} -fPIC")

