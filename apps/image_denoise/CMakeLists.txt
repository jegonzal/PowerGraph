
project(GraphLab)


add_executable(loopybp_denoise loopybp_denoise.cpp)
#TARGET_LINK_LIBRARIES(loopybp_denoise tcmalloc)


#add_executable(gibbs_denoise gibbs_denoise.cpp)
#TARGET_LINK_LIBRARIES(gibbs_denoise tcmalloc)


add_dist_executable(loopybp_denoise_dist loopybp_denoise_dist.cpp)
dist_target_link_libraries(loopybp_denoise_dist tcmalloc)


# add_executable(learning_denoise learning_denoise.cpp)
# TARGET_LINK_LIBRARIES(learning_denoise tcmalloc)





