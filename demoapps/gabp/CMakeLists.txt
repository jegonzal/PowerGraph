project(GraphLab)
add_graphlab_executable(gabp main.cpp advanced_config.cpp unittest.cpp
../../libs/matrixmarket/mmio.cpp logreg.cpp)

