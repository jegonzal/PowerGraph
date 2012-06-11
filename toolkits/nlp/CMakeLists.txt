project(GraphLab)

# Primary executable
add_graphlab_executable(lda_sequential_cgs lda_sequential_cgs.cpp)
add_graphlab_executable(cgs_lda cgs_lda.cpp cgs_lda_vertex_program.cpp)
add_graphlab_executable(cgs_lda2 cgs_lda.cpp cgs_lda_vertex_program2.cpp)








