project(GraphLab)


# search for the python library and headers
find_library(PYTHON27_LIBRARY NAMES python2.7)
find_path(PYTHON_PATH NAMES Python.h PATH_SUFFIXES python python2.7)
find_path(SWIG_BIN NAMES swig)

if(PYTHON27_LIBRARY AND PYTHON_PATH AND SWIG_BIN)
  message(STATUS "Python Dev Environment Found")
  message(STATUS "    Python 2.7 Libary: ${PYTHON27_LIBRARY}")
  message(STATUS "    Python Headers: ${PYTHON_PATH}")
  message(STATUS "    SWIG: ${SWIG_BIN}")
  add_graphlab_executable(py_graphlab_exec py_graphlab.cpp)
  
  include_directories(${PYTHON_PATH})
  target_link_libraries(py_graphlab_exec ${PYTHON27_LIBRARY})

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab_wrap.cxx
    COMMAND ${SWIG_BIN}/swig -c++ -python ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab.i
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/py_graphlab.i
    )
  add_library(py_graphlab SHARED py_graphlab.cpp py_graphlab_wrap.cxx)
  set_target_properties(py_graphlab PROPERTIES PREFIX "_" SUFFIX ".so" LINK_FLAGS -mstackrealign COMPILE_DEFINITIONS PYSHARED_LIB)
  target_link_libraries(py_graphlab ${PYTHON27_LIBRARY} graphlab)
  add_dependencies(py_graphlab graphlab)
  requires_core_deps(py_graphlab)

  copy_file(wrappers.py)
  copy_file(simple_pagerank.py)
  copy_file(als.py)
else()
  message(STATUS "Python Dev Environment Incomplete")
  message(STATUS "    Python 2.7 Libary: ${PYTHON27_LIBRARY}")
  message(STATUS "    Python Headers: ${PYTHON_PATH}")
  message(STATUS "    SWIG: ${SWIG_BIN}")
endif()
