project(netflix++)

add_graphlab_executable(netflix_main netflix_main.cpp)
requires_eigen(netflix_main) # build and attach eigen
