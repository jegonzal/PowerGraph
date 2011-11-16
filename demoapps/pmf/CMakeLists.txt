project(GraphLab)
# include(CheckCXXSourceCompiles)





if(HAS_EIGEN)
  add_graphlab_executable(als als.cpp matrixmarket/mmio.cpp)
endif()
  

  #   add_graphlab_executable(pmf pmf.cpp advanced_config.cpp  ../../libs/matrixmarket/mmio.cpp main.cpp)
  # #  add_library(graphlabpmf STATIC pmf.cpp advanced_config.cpp ../../libs/matrixmarket/mmio.cpp ) 
  # #  ADD_DEFINITIONS(-fPIC)
  # #  add_library(graphlabpmf_pic STATIC pmf.cpp advanced_config.cpp ../../libs/matrixmarket/mmio.cpp )
  #   # test files
  #     add_graphlab_executable(itdiff itdiff.cpp)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest
  #     ${CMAKE_CURRENT_BINARY_DIR}/smalltest COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smallteste
  #     ${CMAKE_CURRENT_BINARY_DIR}/smallteste COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest.out
  #     ${CMAKE_CURRENT_BINARY_DIR}/smalltest.out COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest_eigen.out
  #     ${CMAKE_CURRENT_BINARY_DIR}/smalltest_eigen.out COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/als
  #     ${CMAKE_CURRENT_BINARY_DIR}/als COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/alse
  #     ${CMAKE_CURRENT_BINARY_DIR}/alse COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/wals
  #     ${CMAKE_CURRENT_BINARY_DIR}/wals COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/walse
  #     ${CMAKE_CURRENT_BINARY_DIR}/walse COPYONLY)
  #     configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lanczos2
  #     ${CMAKE_CURRENT_BINARY_DIR}/lanczos2 COPYONLY)

