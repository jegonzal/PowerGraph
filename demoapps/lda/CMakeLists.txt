
project(lda)

# Primary executable
add_executable(lda lda.cpp collapsed_gibbs.cpp corpus.cpp random.cpp)
target_link_libraries(lda boost_program_options itpp)
