project(GraphLab)

# link_libraries(${Boost_LIBRARIES})


subdirs(
#  aggregation
  engine
  graph
  options
  parallel
  scheduler
  util
  serialization
  logger
  ui
#  jni
  )



#build the graphlab library
add_library(graphlab STATIC
  options/command_line_options.cpp
  options/options_map.cpp
  util/timer.cpp
  util/generics/any.cpp
  util/hdfs.cpp
  logger/logger.cpp
  logger/assertions.cpp
  parallel/pthread_tools.cpp
  parallel/thread_pool.cpp
  util/random.cpp
  scheduler/scheduler_list.cpp
  util/net_util.cpp
  util/safe_circular_char_buffer.cpp
  util/fs_util.cpp
  util/memory_info.cpp
  util/tracepoint.cpp
  util/mpi_tools.cpp
  util/web_util.cpp
  rpc/dc_tcp_comm.cpp
  rpc/circular_char_buffer.cpp
  rpc/dc_stream_receive.cpp
  rpc/dc_buffered_stream_send2.cpp
  rpc/dc.cpp
  rpc/reply_increment_counter.cpp
  rpc/dc_init_from_env.cpp
  rpc/dc_init_from_mpi.cpp
  rpc/async_consensus.cpp
  rpc/distributed_event_log.cpp
  rpc/delta_dht.cpp
  ui/mongoose/mongoose.c
  )
requires_core_deps(graphlab)


# if(Sctp-FOUND)
#   set_property(TARGET graphlab APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
#   target_link_libraries(graphlab sctp)
# endif()



#build PIC versions of everything

get_target_property(graphlab_sources graphlab SOURCES)
get_property(graphlab_flags TARGET graphlab PROPERTY COMPILE_FLAGS)
add_library(graphlab_pic STATIC ${graphlab_sources})
requires_core_deps(graphlab_pic)
target_link_libraries(graphlab_pic ${Boost_SHARED_LIBRARIES})


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



