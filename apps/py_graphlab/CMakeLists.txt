project(GraphLab)

add_graphlab_executable(py_graphlab_exec py_graphlab.cpp)
target_link_libraries(py_graphlab_exec python2.7)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab_wrap.cxx
  COMMAND swig -c++ -python ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab.i
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab.i
)
add_library(py_graphlab SHARED py_graphlab.cpp py_graphlab_wrap.cxx)
set_target_properties(py_graphlab PROPERTIES PREFIX "_" SUFFIX ".so" LINK_FLAGS -mstackrealign COMPILE_DEFINITIONS PYSHARED_LIB)
target_link_libraries(py_graphlab python2.7 graphlab)
add_dependencies(py_graphlab graphlab)
requires_core_deps(py_graphlab)

copy_file(wrappers.py)
copy_file(simple_pagerank.py)

