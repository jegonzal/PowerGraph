project(netflix_demo)

add_graphlab_executable(gl3als_netflix_demo gl3als.cpp)
requires_eigen(gl3als_netflix_demo) # build and attach eigen

add_graphlab_executable(sgd_netflix_demo sgd.cpp)
requires_eigen(sgd_netflix_demo) # build and attach eigen

add_graphlab_executable(als_netflix_demo als.cpp)
requires_eigen(als_netflix_demo) # build and attach eigen


add_graphlab_executable(gl3sgd_netflix_demo gl3sgd.cpp)
requires_eigen(gl3sgd_netflix_demo) # build and attach eigen
