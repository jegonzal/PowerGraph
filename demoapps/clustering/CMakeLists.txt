project(GraphLab)

if (ITPP-FOUND)
  add_graphlab_executable(glcluster clustering.cpp advanced_config.cpp read_matrix_market.cpp ../../libs/matrixmarket/mmio.cpp distance.cpp sparse_math.cpp main.cpp kmeans_plus_plus.cpp)
else()
  message(STATUS "ITPP not found. PMF will not be built")
endif()
