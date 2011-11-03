project(GraphLab)

add_graphlab_executable(lda lda.cpp corpus.cpp)

add_graphlab_executable(diff_lda lda.cpp corpus.cpp)
set_target_properties(diff_lda PROPERTIES COMPILE_FLAGS -DDIFFABLE)

add_graphlab_executable(factorized_lda lda.cpp corpus.cpp)
set_target_properties(factorized_lda PROPERTIES COMPILE_FLAGS -DFACTORIZED)



add_graphlab_executable(diff_ss_lda lda.cpp corpus.cpp)
set_target_properties(diff_ss_lda PROPERTIES COMPILE_FLAGS "-DDIFFABLE -DSHAREDSUM")



add_dist2_executable(dht_lda dht_lda.cpp corpus.cpp)

