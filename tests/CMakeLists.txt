project(tests)



link_libraries(${Boost_LIBRARIES})
link_libraries(${GraphLab_LIBRARIES})
include_directories( . )


ADD_CXXTEST(graph_test.cxx)
ADD_CXXTEST(sharedsum_test.cxx)
#ADD_CXXTEST(disk_graph_test.cxx)

add_dist2_executable(delta_dht_test delta_dht_test.cpp)


ADD_CXXTEST(random_test.cxx)
ADD_CXXTEST(factor_test.cxx)
ADD_CXXTEST(small_map_test.cxx)
ADD_CXXTEST(small_set_test.cxx)
ADD_CXXTEST(fast_multinomial_test.cxx)


# ADD_CXXTEST(graphlab_test.cxx)
# ADD_CXXTEST(serializetests.cxx)
# ADD_CXXTEST(thread_tools.cxx)
add_executable(anytests anytests.cpp)
# add_executable(anytests_loader anytests_loader.cpp)

# if (MPI_FOUND)
# add_executable(dc_consensus_test dc_consensus_test.cpp)

# add_executable(rpc_example1 rpc_example1.cpp)
# add_executable(rpc_example2 rpc_example2.cpp)
# add_executable(rpc_example3 rpc_example3.cpp)
# add_executable(rpc_example4 rpc_example4.cpp)
# add_executable(rpc_example5 rpc_example5.cpp)
# add_executable(rpc_example6 rpc_example6.cpp)
# add_executable(rpc_example7 rpc_example7.cpp)

# add_dist2_executable(distributed_dg_construction_test distributed_dg_construction_test.cpp)
# add_dist2_executable(distributed_graph_test distributed_graph_test.cpp)
# endif()

# ADD_CXXTEST(md5test.cxx)

# if(EXPERIMENTAL)
# ADD_CXXTEST(mmap_allocator_test.cxx)
# endif()

# copyfile(runtests.sh)
