project(GraphLab)



# Not supported anymore (single category coem)
# add_executable(coem coem.cpp coemapp.cpp)
# target_link_libraries(coem tcmalloc)


#add_executable(multicoem  multicoemapp.cpp)
#target_link_libraries(multicoem tcmalloc)


# add_executable(multicoemlarge multicoemlarge.cpp multicoemapp.cpp)
# target_link_libraries(multicoemlarge tcmalloc)

add_dist_executable(distmulticoem  distmulticoemapp.cpp)
dist_target_link_libraries(distmulticoem tcmalloc)

add_dist_executable(bandwidth_regression  bandwidth_regression.cpp)
dist_target_link_libraries(bandwidth_regression tcmalloc)


# add_dist_executable(distmulticoem_pull distmulticoem_pull.cpp distmulticoemapp_pull.cpp)
#dist_target_link_libraries(distmulticoem_pull tcmalloc)
