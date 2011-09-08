project(GraphLab)

add_graphlab_executable(lda lda.cpp corpus.cpp)

add_graphlab_executable(diff_lda lda.cpp corpus.cpp)
set_target_properties(diff_lda PROPERTIES COMPILE_FLAGS -DDIFFABLE)



add_dist2_executable(dht_lda dht_lda.cpp corpus.cpp)

