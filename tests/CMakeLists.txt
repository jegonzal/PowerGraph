project(tests)

subdirs(data)

#ADD_CXXTEST(graph_test.cxx)
#ADD_CXXTEST(graph2_test.cxx)



ADD_CXXTEST(random_test.cxx)
# move into toolkit
#ADD_CXXTEST(factor_test.cxx)
ADD_CXXTEST(small_map_test.cxx)
ADD_CXXTEST(small_set_test.cxx)
ADD_CXXTEST(atomic_add_vector.cxx)

ADD_CXXTEST(dense_bitset_test.cxx)

ADD_CXXTEST(serializetests.cxx)
ADD_CXXTEST(thread_tools.cxx)

ADD_CXXTEST(test_lock_free_pool.cxx)
# ADD_CXXTEST(engine_terminator_bench.cxx)
# ADD_CXXTEST(chandy_misra.cxx)
ADD_CXXTEST(lock_free_pushback.cxx)
ADD_CXXTEST(union_find_test.cxx)

ADD_CXXTEST(local_graph_test.cxx)
ADD_CXXTEST(empty_test.cxx)
ADD_CXXTEST(scheduler_test.cxx)

add_graphlab_executable(distributed_graph_test distributed_graph_test.cpp)

add_graphlab_executable(cuckootest cuckootest.cpp)
add_graphlab_executable(dc_consensus_test dc_consensus_test.cpp)
#add_graphlab_executable(distributed_chandy_misra_test distributed_chandy_misra_test.cpp)
add_graphlab_executable(dc_test_sequentialization dc_test_sequentialization.cpp)
add_graphlab_executable(hdfs_test hdfs_test.cpp)
add_graphlab_executable(test_parsers test_parsers.cpp)


add_graphlab_executable(synchronous_engine_test synchronous_engine_test.cpp)
add_graphlab_executable(async_consistent_test async_consistent_test.cpp)

add_graphlab_executable(sfinae_function_test sfinae_function_test.cpp)

# copyfile(runtests.sh)
