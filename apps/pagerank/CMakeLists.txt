project(GraphLab)


add_executable(pagerank pagerank.cpp)
#target_link_libraries(pagerank tcmalloc)

add_jni_library(pagerankjni pagerankjni.cpp)
