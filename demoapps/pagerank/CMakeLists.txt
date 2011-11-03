project(GraphLab)



add_graphlab_executable(pagerank pagerank.cpp graph_loader.cpp)

add_graphlab_executable(factorized_pagerank pagerank.cpp graph_loader.cpp)
set_target_properties(factorized_pagerank PROPERTIES COMPILE_FLAGS -DFACTORIZED)

add_graphlab_executable(diffable_pagerank pagerank.cpp graph_loader.cpp)
set_target_properties(diffable_pagerank PROPERTIES COMPILE_FLAGS -DDIFFABLE)
