project(GraphLab)

# link_libraries(${Boost_LIBRARIES})


subdirs(
 aggregation
 engine
# distributed2
 graph
 gpu
 monitoring
 options
 parallel
 scheduler
 scope
 shared_data
 update_functor
 util
 serialization
 logger
 metrics)

if(YRL_EXPERIMENTAL)
  subdirs(yrl)
endif()

include(metisextern.cmake)

#build graphlab_zoltan if available
## The current distributed build no longer uses zoltan
# IF(ZOLTAN_FOUND)
#   add_library(graphlab_zoltan
#               distributed2/graph/partitioning/adjacency_list.cpp
#               distributed2/graph/partitioning/partitioning_tools.cpp)
# #   TARGET_LINK_LIBRARIES(graphlab_zoltan 
# #     ${MPI_LIBRARY} 
# #     ${MPI_EXTRA_LIBRARY}
# #     ${ZOLTAN_LIBS})
# ENDIF()



#build the graphlab library
if(Sctp-FOUND)
  set(sctp_source rpc/dc_sctp_comm.cpp)
endif()

if (ZLIB_FOUND)
  set(rpc_compressed rpc/dc_stream_receive_z.cpp
                      rpc/dc_buffered_stream_send_expqueue_z.cpp)
endif()

if (MPI_FOUND)
  set(util_mpi_tools util/mpi_tools.cpp)
endif()

IF (EXPERIMENTAL)
  set(mmap_allocator util/mmap_allocator.cpp)
endif()

add_library(graphlab STATIC
  ${metis_sources}
  options/command_line_options.cpp
  options/options_map.cpp
  util/timer.cpp
  util/generics/any.cpp
  graph/disk_atom.cpp
  graph/atom_index_file.cpp
  logger/logger.cpp
  logger/assertions.cpp
  parallel/pthread_tools.cpp
  parallel/thread_pool.cpp
  util/random.cpp
  scheduler/scheduler_list.cpp
  metrics/metrics.cpp
  util/net_util.cpp
  util/safe_circular_char_buffer.cpp
  util/fs_util.cpp
  util/md5.cpp
  ${mmap_allocator}
  factors/discrete_variable.cpp
  factors/binary_factor.cpp
  factors/unary_factor.cpp
  ${util_mpi_tools}
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
  distributed2/distributed_scheduler_list.cpp
  ${rpc_compressed}
)

target_link_libraries(graphlab ${Boost_LIBRARIES})

if (MPI_FOUND)
target_link_libraries(graphlab ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
endif()


if(Sctp-FOUND)
  set_property(TARGET graphlab APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
  target_link_libraries(graphlab sctp)
endif()



# #build distributed graphlab
# add_library(graphlab_dist2 STATIC
#     distributed2/distributed_glshared_base.cpp
#     distributed2/distributed_glshared_manager.cpp
#     distributed2/graph/dgraph_edge_list.cpp)


#build PIC versions of everything

get_target_property(graphlab_sources graphlab SOURCES)
get_property(graphlab_flags TARGET graphlab PROPERTY COMPILE_FLAGS)
add_library(graphlab_pic STATIC ${graphlab_sources})
set_target_properties(graphlab_pic PROPERTIES 
                      COMPILE_FLAGS "${graphlab_flags} -fPIC")

INSTALL(TARGETS 
  graphlab 
  graphlab_pic ARCHIVE DESTINATION lib)

INSTALL(DIRECTORY matlab/  
  DESTINATION include/graphlab/matlab
  FILES_MATCHING PATTERN "*" 
  PATTERN "*.sh" EXCLUDE
)

INSTALL(FILES matlab/get_mex_params.sh
        DESTINATION include/graphlab/matlab
        PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE)


target_link_libraries(graphlab_pic ${Boost_SHARED_LIBRARIES})
