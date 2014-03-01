project(GraphLab)
# include(CheckCXXSourceCompiles)



# Build als
add_graphlab_executable(jacobi jacobi.cpp)
requires_eigen(jacobi) # build and attach eigen


