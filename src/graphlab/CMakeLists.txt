project(GraphLab)

# link_libraries(${Boost_LIBRARIES})


subdirs(
  engine
  graph
  monitoring
  options
  parallel
  scheduler
  context
  update_functor
  util
  serialization
  logger
  metrics
  jni
  )

# Defines METIS_SOURCES
include(extern/metisextern.cmake)


#build the graphlab library
if(Sctp-FOUND)
  set(sctp_source rpc/dc_sctp_comm.cpp)
endif()

if (MPI_FOUND)
  set(util_mpi_tools util/mpi_tools.cpp)
endif()

IF (EXPERIMENTAL)
  set(mmap_allocator util/mmap_allocator.cpp)
endif()

add_library(graphlab STATIC
  ${METIS_SOURCES}
  options/command_line_options.cpp
  options/options_map.cpp
  util/timer.cpp
  util/generics/any.cpp
  util/generics/any_vector.cpp
  graph/disk_atom.cpp
  graph/write_only_disk_atom.cpp
  graph/atom_index_file.cpp
  graph/memory_atom.cpp
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
  util/tracepoint.cpp
  util/event_log.cpp
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
  rpc/dc_buffered_stream_receive.cpp
  rpc/dc.cpp
  rpc/reply_increment_counter.cpp
  rpc/dc_comm_services.cpp
  rpc/dc_init_from_env.cpp
  rpc/dc_init_from_mpi.cpp
  rpc/async_consensus.cpp
  rpc/distributed_event_log.cpp
  rpc/delta_dht.cpp
  )

target_link_libraries(graphlab ${Boost_LIBRARIES})

if (MPI_FOUND)
  target_link_libraries(graphlab ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
endif()


if(Sctp-FOUND)
  set_property(TARGET graphlab APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
  target_link_libraries(graphlab sctp)
endif()



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
