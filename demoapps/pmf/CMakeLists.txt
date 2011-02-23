project(GraphLab)

CHECK_LIBRARY_EXISTS("itpp" "" "" ITPP-FOUND)
if (ITPP-FOUND)
  add_executable(pmf pmf.cpp)
  target_link_libraries(pmf itpp)
endif()
