project(GraphLab)
# include(CheckCXXSourceCompiles)


if (HAS_EIGEN)
  include_directories(../../deps/eigen-eigen-3.0.2/)
  message(status "pointing pmf to eigen directories")
endif()

if (HAS_EIGEN)
  set(lapack_found 1)
elseif(HAS_ITPP)
  set(lapack_found 1)
endif()

if (lapack_found)
  add_graphlab_executable(pmf pmf.cpp advanced_config.cpp  ../../libs/matrixmarket/mmio.cpp main.cpp)
  add_library(graphlabpmf STATIC pmf.cpp advanced_config.cpp ../../libs/matrixmarket/mmio.cpp ) 
  ADD_DEFINITIONS(-fPIC)
  add_library(graphlabpmf_pic STATIC pmf.cpp advanced_config.cpp ../../libs/matrixmarket/mmio.cpp )
  # test files
    add_graphlab_executable(itdiff itdiff.cpp)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest
    ${CMAKE_CURRENT_BINARY_DIR}/smalltest COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smallteste
    ${CMAKE_CURRENT_BINARY_DIR}/smallteste COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest.out
    ${CMAKE_CURRENT_BINARY_DIR}/smalltest.out COPYONLY)
endif()
