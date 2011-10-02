project(GraphLab)
if (HAS_EIGEN)
include_directories(../../deps/eigen-eigen-3.0.2/)
add_definitions(-DHAS_EIGEN)
else()

endif()

add_graphlab_executable(gabp main.cpp advanced_config.cpp unittest.cpp
lasso.cpp
../../libs/matrixmarket/mmio.cpp logreg.cpp)

