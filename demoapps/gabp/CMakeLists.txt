project(GraphLab)
if (HAS_EIGEN)
include_directories(../../deps/eigen-eigen-3.0.2/)
add_definitions(-DHAS_EIGEN)
else()

endif()

if (HAS_EIGEN)
  set( lapack_found  1)
elseif (ITPP-FOUND )
  set( lapack_found  1)
endif()

if (lapack_found)
add_graphlab_executable(gabp main.cpp advanced_config.cpp unittest.cpp
lasso.cpp gamp.cpp
../../libs/matrixmarket/mmio.cpp logreg.cpp)
endif()
