project(GraphLab)

# Primary executable
add_graphlab_executable(lda_sequential_cgs lda_sequential_cgs.cpp)
add_graphlab_executable(cgs_lda_gas cgs_lda.cpp cgs_lda_vertex_program_gas.cpp)
add_graphlab_executable(cgs_lda cgs_lda.cpp cgs_lda_vertex_program.cpp)








