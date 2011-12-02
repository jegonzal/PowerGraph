project(GraphLab)
add_graphlab_executable(jacobi jacobi.cpp ../shared/mmio.cpp)
add_graphlab_executable(jacobi2 jacobi2.cpp ../shared/mmio.cpp)
add_graphlab_executable(gabp gabp.cpp ../shared/mmio.cpp)

