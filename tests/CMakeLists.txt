project(tests)

link_libraries(${GraphLab_LIBRARIES})

include_directories( . )

ADD_CXXTEST(graphtest.cxx)

ADD_CXXTEST(serializetests.cxx)
ADD_CXXTEST(foreachtest.cxx)
ADD_CXXTEST(logtest.cxx)
ADD_CXXTEST(circqueuetest.cxx)
# ADD_CXXTEST(basictasktests.cxx)
# ADD_CXXTEST(basicgraphlabintegrationtest.cxx)
ADD_CXXTEST(sampletreetest.cxx)
# ADD_CXXTEST(set_scheduler_test.cxx)
ADD_CXXTEST(shared_data_test.cxx)
ADD_CXXTEST(tlstest.cxx)
ADD_CXXTEST(factor_tests.cxx)
# ADD_CXXTEST(fastmultinomialtest.cxx)
# ADD_CXXTEST(edge_selector_test.cxx)
# ADD_CXXTEST(randomtest.cxx)
# ADD_CXXTEST(vertex_set_test.cxx)
ADD_CXXTEST(dense_bitset_test.cxx)
ADD_CXXTEST(timertest.cxx)
# ADD_CXXTEST(executionplantest.cxx)
add_executable(programoptions programoptions.cpp)
# add_executable(test_optimize test_optimize.cpp)
add_dist_executable(distinittest distinittest.cpp)
add_dist_executable(distmsgtest distmsgtest.cpp)
add_executable(anytests anytests.cpp)
add_executable(anytests_loader anytests_loader.cpp)
add_executable(coretest coretest.cpp)
add_dist_executable(clonedgraphtest clonedgraphtest.cpp)
add_dist_executable(dist_lm_test dist_lm_test.cpp)
add_dist_executable(dist_lm_test2 dist_lm_test2.cpp)
add_dist_executable(distdht distdht.cpp)
add_dist_executable(dist_lm_test_block dist_lm_test_block.cpp)
add_dist_executable(dist_engine_test dist_engine_test.cpp)
add_dist_executable(dist_engine_test2 dist_engine_test2.cpp)
add_dist_executable(dist_engine_test3 dist_engine_test3.cpp)
add_dist_executable(dist_term_test dist_term_test.cpp)
add_dist_executable(distributed_shared_data_test distributed_shared_data_test.cpp)

# cuda_add_executable(cuda_test cuda_test.cpp)

# add_executable(itppserialize itppserialize.cpp)
#target_link_libraries(itppserialize  itpp)


add_executable(partransform partransform.cpp)

add_executable(clopts command_line_options.cpp)


# serialization tests 
add_executable(integercompress integercompress.cpp)
add_executable(integercompress_bench integercompress_bench.cpp)
add_executable(podtest podtest.cpp)
add_executable(podtest_loader podtest_loader.cpp)

#add_executable(circular_test circular_test.cpp)
#add_executable(simple_function_launch simple_function_launch.cpp)
#add_executable(dctest dctest.cpp)
#add_executable(dcobjtest dcobjtest.cpp)
#add_executable(is_rpc_call_test is_rpc_call_test.cpp)
#add_executable(dc_barrier_test dc_barrier_test.cpp)
#add_executable(dhttest dhttest.cpp)
#add_executable(dht_performance_test dht_performance_test.cpp)

#add_executable(example1 example1.cpp)
#add_executable(example2 example2.cpp)
#add_executable(example3 example3.cpp)
#add_executable(example4 example4.cpp)
#add_executable(example5 example5.cpp)
#add_executable(example6 example6.cpp)
#add_executable(example7 example7.cpp)
#add_executable(example8 example8.cpp)



