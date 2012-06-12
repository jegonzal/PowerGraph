project(GraphLab)

# Primary executable
add_graphlab_executable(lda_sequential_cgs lda_sequential_cgs.cpp)
add_graphlab_executable(cgs_ldas cgs_lda.cpp cgs_lda_common.cpp)
add_graphlab_executable(fast_cgs_ldas fast_cgs_lda.cpp cgs_lda_common.cpp)









