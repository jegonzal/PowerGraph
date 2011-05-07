project(tests)
link_libraries(${Boost_LIBRARIES})
link_libraries(${GraphLab_LIBRARIES})
include_directories( . )


ADD_CXXTEST(graph_test.cxx)
ADD_CXXTEST(disk_graph_test.cxx)
target_link_libraries(disk_graph_test.cxxtest kyotocabinet)

ADD_CXXTEST(randomtest.cxx)
ADD_CXXTEST(graphlab_test.cxx)
ADD_CXXTEST(serializetests.cxx)

add_executable(anytests anytests.cpp)
add_executable(anytests_loader anytests_loader.cpp)
add_executable(dc_consensus_test dc_consensus_test.cpp)

ADD_CXXTEST(md5test.cxx)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh
    ${CMAKE_CURRENT_BINARY_DIR}/runtests.sh COPYONLY)
