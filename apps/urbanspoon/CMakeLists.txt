project(urbanspoon)

add_graphlab_executable(main main.cpp)
requires_eigen(main) # build and attach eigen

add_graphlab_executable(data_split data_split.cpp)
