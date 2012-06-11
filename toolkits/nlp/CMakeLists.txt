project(GraphLab)

# Primary executable
add_graphlab_executable(lda_sequential_cgs lda_sequential_cgs.cpp)
add_graphlab_executable(cgs_ldas cgs_lda.cpp 
  fast_cgs_lda_vertex_program.cpp cgs_lda_vertex_program.cpp)








