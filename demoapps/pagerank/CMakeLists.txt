project(GraphLab)



add_graphlab_executable(pagerank pagerank.cpp graph_loader.cpp)

add_graphlab_executable(pagerank_factorized pagerank.cpp graph_loader.cpp)

set_target_properties(pagerank_factorized PROPERTIES COMPILE_FLAGS -DFACTORIZED)
