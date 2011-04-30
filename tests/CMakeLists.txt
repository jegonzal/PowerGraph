project(tests)
link_libraries(${Boost_LIBRARIES})
link_libraries(${GraphLab_LIBRARIES})
include_directories( . )


ADD_CXXTEST(graph_test.cxx)
