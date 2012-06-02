project(GraphLab)
# include(CheckCXXSourceCompiles)


if (HAS_EIGEN)
  include_directories(../../deps/eigen-eigen-3.0.2/)
  message(status "pointing pmf to eigen directories")
  add_definitions(-DHAS_EIGEN)
endif()

if (HAS_EIGEN)
  set(lapack_found 1)
elseif(ITPP-FOUND)  
  set(lapack_found 1)
endif()

if (lapack_found)
  add_graphlab_executable(pmf pmf.cpp advanced_config.cpp
../../libs/matrixmarket/mmio.cpp main.cpp)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/als
    ${CMAKE_CURRENT_BINARY_DIR}/als COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/alse
    ${CMAKE_CURRENT_BINARY_DIR}/alse COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/wals
    ${CMAKE_CURRENT_BINARY_DIR}/wals COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/walse
    ${CMAKE_CURRENT_BINARY_DIR}/walse COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lanczos2
    ${CMAKE_CURRENT_BINARY_DIR}/lanczos2 COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lanczos2t
    ${CMAKE_CURRENT_BINARY_DIR}/lanczos2t COPYONLY)
endif()
