project(GraphLab)
# include(CheckCXXSourceCompiles)


if (ITPP-FOUND)
  add_graphlab_executable(pmf pmf.cpp)

  # test files
    add_graphlab_executable(itdiff itdiff.cpp)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest
    ${CMAKE_CURRENT_BINARY_DIR}/smalltest COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smallteste
    ${CMAKE_CURRENT_BINARY_DIR}/smallteste COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/smalltest.out
    ${CMAKE_CURRENT_BINARY_DIR}/smalltest.out COPYONLY)
else()
  message(STATUS "ITPP not found. PMF will not be built")
endif()
