project(Extensions)

if(CPP11)

add_library(graphlab_extension STATIC 
            MurmurHash3.cpp 
            extension.cpp 
            extension_graph.cpp
            extension_pagerank.cpp)
add_dependencies(graphlab_extension graphlab)
target_link_libraries(graphlab_extension graphlab)

add_extension_executable(pagerank_extension_driver
                         pagerank_extension_driver.cpp) 

endif()
