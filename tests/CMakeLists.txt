project(tests)
subdirs(data)

ADD_CXXTEST(random_test.cxx)

# move into toolkit
#ADD_CXXTEST(factor_test.cxx)
ADD_CXXTEST(small_map_test.cxx)
ADD_CXXTEST(small_set_test.cxx)

ADD_CXXTEST(dense_bitset_test.cxx)
ADD_CXXTEST(serializetests.cxx)
ADD_CXXTEST(thread_tools.cxx)

ADD_CXXTEST(test_lock_free_pool.cxx)
ADD_CXXTEST(lock_free_pushback.cxx)
ADD_CXXTEST(union_find_test.cxx)

ADD_CXXTEST(empty_test.cxx)
# ADD_CXXTEST(scheduler_test.cxx)

ADD_CXXTEST(csr_storage_test.cxx)
ADD_CXXTEST(local_graph_test.cxx)
add_graphlab_executable(distributed_graph_test distributed_graph_test.cpp)
add_graphlab_executable(distributed_ingress_test distributed_ingress_test.cpp)

add_graphlab_executable(cuckootest cuckootest.cpp)
add_graphlab_executable(dc_consensus_test dc_consensus_test.cpp)
add_graphlab_executable(distributed_chandy_misra_test distributed_chandy_misra_test.cpp)
add_graphlab_executable(dc_fiber_consensus_test dc_fiber_consensus_test.cpp)
add_graphlab_executable(dc_test_sequentialization dc_test_sequentialization.cpp)
add_graphlab_executable(hdfs_test hdfs_test.cpp)
add_graphlab_executable(test_parsers test_parsers.cpp)

add_graphlab_executable(synchronous_engine_test synchronous_engine_test.cpp)
add_graphlab_executable(async_consistent_test async_consistent_test.cpp)

add_graphlab_executable(sfinae_function_test sfinae_function_test.cpp)

add_test(synchronous_engine_test synchronous_engine_test)
add_test(async_consistent_test async_consistent_test)

# copyfile(runtests.sh)

add_graphlab_executable(mini_web_server mini_web_server.cpp)

add_graphlab_executable(test_vertex_set test_vertex_set.cpp)

add_test(test_vertex_set test_vertex_set)
add_graphlab_executable(arbitrary_signal_test arbitrary_signal_test.cpp)


add_graphlab_executable(sort_test sort_test.cpp)

add_graphlab_executable(hopscotch_test hopscotch_test.cpp)

add_graphlab_executable(fiber_test fiber_test.cpp)
add_graphlab_executable(fibo_fiber_test fibo_fiber_test.cpp)
