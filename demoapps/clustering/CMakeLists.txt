project(GraphLab)
if (HAS_EIGEN)
  include_directories(../../deps/eigen-eigen-3.0.2/)
  add_definitions(-DHAS_EIGEN)
endif()

set(OMP_SUPPORT 0)
  if (OMP_SUPPORT)
    add_definitions(-omp)
    add_definitions(-DOMP_SUPPORT)
  endif()

if (HAS_EIGEN)
  set( lapack_found 1)
elseif(ITPP-FOUND)
  set( lapack_found 1)
endif()

if (lapack_found)
  add_graphlab_executable(glcluster clustering.cpp advanced_config.cpp 
read_matrix_market.cpp ../../libs/matrixmarket/mmio.cpp distance.cpp   
sparse_math.cpp main.cpp kmeans_plus_plus.cpp dumpcluster.cpp tfidf.cpp
lda.cpp gamma.cpp kcores.cpp)
endif()

