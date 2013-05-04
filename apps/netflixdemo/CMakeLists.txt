project(netflix_demo)

# liblbfgs =================================================================
ExternalProject_Add(liblbfgs
  PREFIX ${CMAKE_SOURCE_DIR}/deps/liblbfgs
  URL https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz
  URL_MD5 2a46da6c4014d6b1e8a8790a93edffbb
  INSTALL_DIR ${CMAKE_SOURCE_DIR}/deps/local
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-sse2 
  BUILD_IN_SOURCE 1
  BUILD_COMMAND  ""
  INSTALL_COMMAND make install
  )

macro(requires_lbfgs NAME)
  add_dependencies(${NAME} liblbfgs)
  target_link_libraries(${NAME} lbfgs)
endmacro(requires_lbfgs)


add_graphlab_executable(gl3als_netflix_demo gl3als.cpp)
requires_eigen(gl3als_netflix_demo) # build and attach eigen

add_graphlab_executable(sgd_netflix_demo sgd.cpp)
requires_eigen(sgd_netflix_demo) # build and attach eigen

add_graphlab_executable(als_netflix_demo als.cpp)
requires_eigen(als_netflix_demo) # build and attach eigen


add_graphlab_executable(gl3sgd_netflix_demo gl3sgd.cpp)
requires_eigen(gl3sgd_netflix_demo) # build and attach eigen
