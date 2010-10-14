project(GraphLab)


add_executable(PMF pmf.cpp)
target_link_libraries(PMF tcmalloc)
target_link_libraries(PMF itpp)
#add_executable(PMF2 pmf2.cpp)
#target_link_libraries(PMF2 tcmalloc)
#target_link_libraries(PMF2 itpp_debug)
#add_executable(PMF3 pmf3.cpp)
#target_link_libraries(PMF3 tcmalloc)
#target_link_libraries(PMF3 itpp_debug)
#add_dist_executable(PMF4 pmf4.cpp)
#target_link_libraries(PMF4 tcmalloc)
#target_link_libraries(PMF4 itpp)
#add_dist_executable(PMF5 pmf5.cpp)
##target_link_libraries(PMF5 tcmalloc)
#target_link_libraries(PMF5 itpp)

add_dist_executable(PMF5push pmf5push.cpp)
dist_target_link_libraries(PMF5push tcmalloc)
dist_target_link_libraries(PMF5push itpp)
#target_link_libraries(PMF5push tcmalloc)

add_dist_executable(PMF5pushbest pmf5pushbest.cpp)
dist_target_link_libraries(PMF5pushbest tcmalloc)
dist_target_link_libraries(PMF5pushbest itpp)



#
#add_executable(DPMF dpmf.cpp)
#target_link_libraries(DPMF itpp)
