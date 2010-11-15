project(GraphLab)

add_executable(create_denoise_graph create_denoise_graph.cpp)

# add_executable(sequential_denoise sequential_denoise.cpp)


add_executable(pgibbs_denoise pgibbs_denoise.cpp)
TARGET_LINK_LIBRARIES(pgibbs_denoise tcmalloc)

add_executable(pgibbs_protein pgibbs_protein.cpp)
TARGET_LINK_LIBRARIES(pgibbs_protein tcmalloc)


add_executable(priority_queue_test priority_queue_test.cpp)


add_executable(pgibbs_ising pgibbs_ising.cpp)
