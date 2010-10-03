project(GraphLab)



add_executable(makesuperpix makesuperpix.cpp)
target_link_libraries(makesuperpix jpeg tcmalloc)

add_executable(makesuperpix2 makesuperpix2.cpp)
target_link_libraries(makesuperpix2 jpeg tcmalloc)

add_executable(makesuperpix3 makesuperpix3.cpp)
target_link_libraries(makesuperpix3 jpeg tcmalloc)

add_executable(makesuperpix4 makesuperpix4.cpp)
target_link_libraries(makesuperpix4 jpeg tcmalloc)

add_executable(cosegment cosegment.cpp)
target_link_libraries(cosegment jpeg)


add_executable(makefeatures makefeatures.cpp)
target_link_libraries(makefeatures jpeg)

add_dist_executable(cosegment_dist cosegment_dist.cpp)
dist_target_link_libraries(cosegment_dist jpeg tcmalloc)

add_dist_executable(cosegment_dist2 cosegment_dist2.cpp)
dist_target_link_libraries(cosegment_dist2 jpeg tcmalloc)

add_dist_executable(cosegment_dist3 cosegment_dist3.cpp)
dist_target_link_libraries(cosegment_dist3 jpeg tcmalloc)

add_executable(combinebins combinebins.cpp)
add_executable(cosegment3 cosegment3.cpp)
target_link_libraries(cosegment3 jpeg tcmalloc)


add_dist_executable(makesuperpix4_dist makesuperpix4_dist.cpp)
dist_target_link_libraries(makesuperpix4_dist jpeg tcmalloc)