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
  logger/backtrace.cpp
  parallel/pthread_tools.cpp
  # parallel/qthread_tools.cpp
  parallel/thread_pool.cpp
  parallel/fiber_control.cpp
  parallel/fiber_group.cpp
  util/random.cpp
  scheduler/scheduler_list.cpp
  scheduler/fifo_scheduler.cpp
  scheduler/priority_scheduler.cpp
  scheduler/sweep_scheduler.cpp
  scheduler/queued_fifo_scheduler.cpp
  util/net_util.cpp
  util/safe_circular_char_buffer.cpp
  util/fs_util.cpp
  util/memory_info.cpp
  util/tracepoint.cpp
  util/mpi_tools.cpp
  util/web_util.cpp
  util/inplace_lf_queue.cpp
  zookeeper/zookeeper_common.cpp
  zookeeper/key_value.cpp
  zookeeper/server_list.cpp
  rpc/dc_tcp_comm.cpp
  rpc/circular_char_buffer.cpp
  rpc/dc_stream_receive.cpp
  rpc/dc_buffered_stream_send2.cpp
  rpc/dc.cpp
  rpc/request_reply_handler.cpp
  rpc/dc_init_from_env.cpp
  rpc/dc_init_from_mpi.cpp
  rpc/dc_init_from_zookeeper.cpp
  rpc/async_consensus.cpp
  rpc/fiber_async_consensus.cpp
  rpc/distributed_event_log.cpp
  rpc/delta_dht.cpp
  rpc/thread_local_send_buffer.cpp
  ui/mongoose/mongoose.cpp
  ui/metrics_server.cpp
  rpc/get_current_process_hash.cpp
  )
requires_core_deps(graphlab)


# if(Sctp-FOUND)
#   set_property(TARGET graphlab APPEND PROPERTY COMPILE_FLAGS -DHAS_SCTP)
#   target_link_libraries(graphlab sctp)
# endif()




INSTALL(TARGETS 
  graphlab ARCHIVE DESTINATION lib)




