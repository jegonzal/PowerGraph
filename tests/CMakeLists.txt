project(tests)

subdirs(data)

ADD_CXXTEST(graph_test.cxx)
ADD_CXXTEST(graph2_test.cxx)



ADD_CXXTEST(random_test.cxx)
ADD_CXXTEST(factor_test.cxx)
ADD_CXXTEST(small_map_test.cxx)
ADD_CXXTEST(small_set_test.cxx)
ADD_CXXTEST(fast_multinomial_test.cxx)


# ADD_CXXTEST(graphlab_test.cxx)
ADD_CXXTEST(serializetests.cxx)
ADD_CXXTEST(thread_tools.cxx)

ADD_CXXTEST(test_lock_free_pool.cxx)
ADD_CXXTEST(engine_terminator_bench.cxx)
ADD_CXXTEST(md5test.cxx)
ADD_CXXTEST(lock_free_queues.cxx)

add_graphlab_executable(integer_serialize_validate integer_serialize_validate.cpp)
add_graphlab_executable(integer_serialize_bench integer_serialize_bench.cpp)
add_graphlab_executable(test_graph_ops test_graph_ops.cpp)
add_graphlab_executable(dc_consensus_test dc_consensus_test.cpp)
add_graphlab_executable(delta_dht_test delta_dht_test.cpp)



# copyfile(runtests.sh)
