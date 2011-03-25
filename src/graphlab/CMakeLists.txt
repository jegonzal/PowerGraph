project(GraphLab)

# link_libraries(${Boost_LIBRARIES})


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
 logger
 metrics
 jni)


if (MPI_FOUND)
  set(util_mpi_tools util/mpi_tools.cpp)
endif()

add_library(graphlab_util STATIC
  util/timer.cpp
  util/generics/any.cpp
  util/command_line_options.cpp
  logger/logger.cpp
  logger/assertions.cpp
  parallel/pthread_tools.cpp
  schedulers/scheduler_list.cpp
  metrics/metrics.cpp
  util/net_util.cpp
  util/safe_circular_char_buffer.cpp
  util/fs_util.cpp
  ${util_mpi_tools}
)

if (MPI_FOUND)
  target_link_libraries(graphlab_util ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
endif()



if(Sctp-FOUND)
  set(sctp_source rpc/dc_sctp_comm.cpp)
endif()

if (ZLIB_FOUND)
  set(rpc_compressed rpc/dc_stream_receive_z.cpp
                      rpc/dc_buffered_stream_send_expqueue_z.cpp)
endif()
  
add_library(graphlab_rpc STATIC
  rpc/dc_comm_base.cpp
  rpc/dc_tcp_comm.cpp
  ${sctp_source}
  rpc/circular_char_buffer.cpp
  rpc/dc_stream_send.cpp
  rpc/dc_stream_receive.cpp
  rpc/dc_buffered_stream_send.cpp
  rpc/dc_buffered_stream_send_expqueue.cpp
  rpc/dc_buffered_stream_send_expqueue2.cpp
  rpc/dc_buffered_stream_receive.cpp
  rpc/dc.cpp
  rpc/reply_increment_counter.cpp
  rpc/dc_comm_services.cpp
  rpc/dc_init_from_env.cpp
  rpc/dc_init_from_mpi.cpp
  rpc/async_consensus.cpp
  ${rpc_compressed}
  rpc/dc_buffered_stream_send_expqueue_lz.cpp
  rpc/dc_stream_receive_lz.cpp)


IF(ZOLTAN_FOUND)
  add_library(graphlab_zoltan
              distributed2/graph/partitioning/adjacency_list.cpp
              distributed2/graph/partitioning/partitioning_tools.cpp)
#   TARGET_LINK_LIBRARIES(graphlab_zoltan 
#     ${MPI_LIBRARY} 
#     ${MPI_EXTRA_LIBRARY}
#     ${ZOLTAN_LIBS})
ENDIF()


if(Sctp-FOUND)
  set_property(TARGET graphlab_rpc APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
endif()

target_link_libraries(graphlab_rpc graphlab_util quicklz)


add_library(graphlab STATIC
  factors/table_factor.cpp
  factors/binary_factor.cpp
  factors/unary_factor.cpp
  schedulers/scheduler_options.cpp)
  

if(Sctp-FOUND)
  target_link_libraries(graphlab_rpc sctp)
endif()

target_link_libraries(graphlab graphlab_metis graphlab_rpc graphlab_util)



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
  ${MPI_EXTRA_LIBRARY})

   TARGET_LINK_LIBRARIES(graphlab_dist ${MPI_LIBRARY}
     ${MPI_EXTRA_LIBRARY} )
   set_property(TARGET graphlab_dist 
     APPEND PROPERTY COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
 endif(MPI_FOUND)


add_library(graphlab_dist2 STATIC
    distributed2/graph/partitioning/adjacency_list.cpp      
    distributed2/graph/atom_index_file.cpp
    distributed2/distributed_glshared_base.cpp
    distributed2/distributed_glshared_manager.cpp)



target_link_libraries(graphlab ${Boost_LIBRARIES})
target_link_libraries(graphlab_util ${Boost_LIBRARIES})



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

INSTALL(TARGETS 
  graphlab 
  graphlab_util 
  graphlab_pic 
  graphlab_util_pic ARCHIVE DESTINATION lib)

INSTALL(DIRECTORY matlab/  
  DESTINATION include/graphlab/matlab
  FILES_MATCHING PATTERN "*" 
  PATTERN "*.sh" EXCLUDE
)

INSTALL(FILES matlab/get_mex_params.sh
        DESTINATION include/graphlab/matlab
        PERMISSIONS OWNER_WRITE WORLD_READ WORLD_EXECUTE)


target_link_libraries(graphlab_pic ${Boost_SHARED_LIBRARIES})
target_link_libraries(graphlab_util_pic ${Boost_SHARED_LIBRARIES})
