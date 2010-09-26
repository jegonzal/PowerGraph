project(GraphLab)



# add_executable(sequential_denoise sequential_denoise.cpp)


add_executable(alchemy_to_clique_graph alchemy_to_clique_graph.cpp)

add_executable(pgibbs_factor_tests pgibbs_factor_tests.cpp)

add_executable(pgibbs_factor_denoise pgibbs_factor_denoise.cpp)
TARGET_LINK_LIBRARIES(pgibbs_factor_denoise tcmalloc)

add_executable(pgibbs_factor_alchemy pgibbs_factor_alchemy.cpp)
TARGET_LINK_LIBRARIES(pgibbs_factor_alchemy tcmalloc)

# add_executable(pgibbs_protein pgibbs_protein.cpp)
# TARGET_LINK_LIBRARIES(pgibbs_protein tcmalloc)


# add_executable(priority_queue_test priority_queue_test.cpp)







