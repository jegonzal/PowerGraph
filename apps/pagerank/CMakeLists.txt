project(GraphLab)


add_executable(pagerank pagerank.cpp)
#add_executable(tsv_to_graphlab_bin tsv_to_graphlab_bin.cpp)
target_link_libraries(pagerank tcmalloc)




#add_mex_library(pagerankmex pagerankapp.cpp  pagerankmex.cpp)
