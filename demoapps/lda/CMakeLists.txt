project(GraphLab)

add_graphlab_executable(lda lda.cpp corpus.cpp)

add_dist2_executable(dht_lda dht_lda.cpp corpus.cpp)

