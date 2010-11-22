project(GraphLab)



add_executable(svm svm.cpp)
target_link_libraries(svm tcmalloc)
