project(GraphLab)

# We require the most recent version of cmake and automatically
# install the correct version when using the cmake lists
cmake_minimum_required(VERSION 2.8)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
# Preprocessor definition values are now escaped automatically.
cmake_policy(SET CMP0005 NEW)

# Determine where additional GraphLab specific cmake modules are
# defined
set(CMAKE_MODULE_PATH ${GraphLab_SOURCE_DIR}/cmake)

# TODO: what is this used by
set(BUILD_SHARED_LIBS ON)

include(CheckLibraryExists)
include(CheckFunctionExists)
include(ExternalProject)
include(CheckCXXSourceCompiles)

# check we are on a 64 bit machine (else fail)
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
   message( SEND_ERROR "GraphLab does not support 32 bit systems. Please switch to 64 bit system and try again" )
endif ( CMAKE_SIZEOF_VOID_P EQUAL 4 )

# # Display information about cmake variables
# include(CMakePrintSystemInformation)

# set include path
include_directories(
  ${GraphLab_SOURCE_DIR}/src
  ${GraphLab_SOURCE_DIR}/cxxtest
  ${GraphLab_SOURCE_DIR}/deps/local/include)

# set link path
link_directories(${GraphLab_SOURCE_DIR}/deps/local/lib)

add_definitions(-DUSE_DYNAMIC_LOCAL_GRAPH)

if(NO_OPENMP)
  set(OPENMP_C_FLAGS "")
  set(OPENMP_LIBRARIES "")
  add_definitions(-D__NO_OPENMP__)
else()
  set(OPENMP_C_FLAGS "-fopenmp")
  set(OPENMP_LIBRARIES "gomp")
endif()

if (APPLE)
link_libraries(pthread ${OPENMP_LIBRARIES})
else()
link_libraries(pthread rt ${OPENMP_LIBRARIES})
endif()

set(CPP11_FLAGS "")
if(CPP11)
  message(STATUS "C++11 Enabled")
  set(CPP11_FLAGS "-std=c++11 -Wno-enum-compare -Wno-conversion-null")
endif()

if(VID32)
  message(STATUS "Using 32bit vertex id types")
  add_definitions(-DUSE_VID32)
endif()


# Shared compiler flags used by all builds (debug, profile, release)
set(COMPILER_FLAGS "-Wall -g ${CPP11_FLAGS} ${OPENMP_C_FLAGS}" CACHE STRING "common compiler options")

# Disable address space randomization for OSX lion and above
if (APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Xlinker -no_pie")
endif()

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,${GraphLab_SOURCE_DIR}/deps/local/lib")

# Set subdirectories
subdirs(src tests demoapps toolkits)
if(EXPERIMENTAL)
  if (IS_DIRECTORY ${GraphLab_SOURCE_DIR}/experimental)
    subdirs(experimental)
  endif()
endif()

if (IS_DIRECTORY ${GraphLab_SOURCE_DIR}/apps)
  subdirs(apps)
endif()

# Set installation rules
install(DIRECTORY src/
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp"  PATTERN "*.h"
  PATTERN ".svn" EXCLUDE
  )
install(CODE "execute_process(COMMAND ./scripts/install_deps.sh
  ${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR})")

# set(EXTERNAL_PROJECT_FLAGS
#   "-I${GraphLab_SOURCE_DIR}/deps/local/include -L${GraphLab_SOURCE_DIR}/deps/local/lib")



### ===================================================================
# ## Download News update
# if(NOT DISABLE_NEWS)
#   message(STATUS
#     "\n"
#     "===============================================================\n"
#     "Downloading the news from graphlab.org: ")
#   file(DOWNLOAD "http://graphlab.org/news.txt"
#     ${GraphLab_SOURCE_DIR}/new_news.txt
#     TIMEOUT 3
#     STATUS download_status)
#   # Track usage v2 news url
#   file(DOWNLOAD "http://c.statcounter.com/7396904/0/a1b29b21/1/"
#     ${GraphLab_SOURCE_DIR}/temporary_download_counter.tmp
#     TIMEOUT 3
#     STATUS dlstatus)
#   file(REMOVE ${GraphLab_SOURCE_DIR}/temporary_download_counter.tmp)
#   file(READ ${GraphLab_SOURCE_DIR}/new_news.txt news_update)
#   message(STATUS
#     "Finished.\n"
#     ${news_update}
#     "\n"
#     "===============================================================")
# endif()






## ============================================================================
## ============================================================================
## ============================================================================
# Build external dependencies

# Build eigen =================================================================
# TODO: Move to toolkits.
# While eigen is not used in the core GraphLab library we found it
# useful to have a matrix library. For now we download eigen from
# our servers since the eigen servers use SSL which is not supported
# by cmake DOWNLOAD_COMMAND hg clone
# https://bitbucket.org/eigen/eigen/
ExternalProject_Add(eigen
  PREFIX ${GraphLab_SOURCE_DIR}/deps/eigen
  URL http://bitbucket.org/eigen/eigen/get/3.1.2.tar.bz2
  URL_MD5 e9c081360dde5e7dcb8eba3c8430fde2
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND cp -r Eigen unsupported <INSTALL_DIR>/
  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local/include)
add_definitions(-DHAS_EIGEN)

macro(requires_eigen NAME)
  add_dependencies(${NAME} eigen)
endmacro(requires_eigen)


# libbz  =================================================================
ExternalProject_Add(libbz2
  PREFIX ${GraphLab_SOURCE_DIR}/deps/libbz2
  URL http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
  URL_MD5 00b516f4704d4a7cb50a1d97e6e8e15b
  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
  CONFIGURE_COMMAND ""
  PATCH_COMMAND patch -N -p0 -i ${GraphLab_SOURCE_DIR}/patches/libbz2_fpic.patch || true
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make install PREFIX=<INSTALL_DIR>
  INSTALL_COMMAND "" )




# Build Boost =================================================================
# GraphLab relies heavily on the boost library for parsing program
# options, filesystem support, random number generation, and
# unordered containers.
ExternalProject_Add(boost
  PREFIX ${GraphLab_SOURCE_DIR}/deps/boost
  URL "http://tcpdiag.dl.sourceforge.net/project/boost/boost/1.53.0/boost_1_53_0.tar.gz" 
  URL_MD5 57a9e2047c0f511c4dfcf00eb5eb2fbb
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND
  ./bootstrap.sh
  --with-libraries=filesystem
  --with-libraries=program_options
  --with-libraries=system
  --with-libraries=iostreams
  --with-libraries=date_time
  --with-libraries=random
  --with-libraries=context
  --prefix=<INSTALL_DIR>
  BUILD_COMMAND
  C_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
  CPLUS_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
  LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib
  ./b2 install link=static variant=release threading=multi runtime-link=static
  INSTALL_COMMAND ""
  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local )
set(BOOST_ROOT ${GraphLab_SOURCE_DIR}/deps/local )
set(BOOST_LIBS_DIR ${GraphLab_SOURCE_DIR}/deps/local/lib)
set(Boost_LIBRARIES
  ${BOOST_LIBS_DIR}/libboost_filesystem.a
  ${BOOST_LIBS_DIR}/libboost_program_options.a
  ${BOOST_LIBS_DIR}/libboost_system.a
  ${BOOST_LIBS_DIR}/libboost_iostreams.a
  ${BOOST_LIBS_DIR}/libboost_context.a
  ${BOOST_LIBS_DIR}/libboost_date_time.a)
add_dependencies(boost libbz2)
message(STATUS "Boost libs: " ${Boost_LIBRARIES})
# add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
# add_definitions(-DBOOST_ALL_DYN_LINK)
# set(Boost_SHARED_LIBRARIES "")
foreach(blib ${Boost_LIBRARIES})
  message(STATUS "Boost libs: " ${blib})
  string(REGEX REPLACE "\\.a$" ${CMAKE_SHARED_LIBRARY_SUFFIX} bout ${blib})
  message(STATUS "Boost dyn libs: " ${bout})
  set(Boost_SHARED_LIBRARIES ${Boost_SHARED_LIBRARIES} ${bout})
endforeach()
message(STATUS "Boost Shared libs: " ${Boost_SHARED_LIBRARIES})

if(NOT NO_TCMALLOC)
    if(APPLE)
    set (tcmalloc_shared "--enable-shared=yes")
    else()
    set (tcmalloc_shared "--enable-shared=no")
    endif()

    # TCMalloc  ===================================================================
    # We use tcmalloc for improved memory allocation performance
    ExternalProject_Add(libtcmalloc
    PREFIX ${GraphLab_SOURCE_DIR}/deps/tcmalloc
    # Some users can't access domain googlecode.com ,This is a spare URL
    # URL http://sourceforge.jp/projects/sfnet_gperftools.mirror/downloads/gperftools-2.0.tar.gz
    URL http://gperftools.googlecode.com/files/gperftools-2.0.tar.gz
    URL_MD5 13f6e8961bc6a26749783137995786b6
    PATCH_COMMAND patch -N -p0 -i ${GraphLab_SOURCE_DIR}/patches/tcmalloc.patch || true
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --enable-frame-pointers --prefix=<INSTALL_DIR> ${tcmalloc_shared}
    INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)
    #link_libraries(tcmalloc)
    set(TCMALLOC-FOUND 1)
    add_definitions(-DHAS_TCMALLOC)
endif()



# Lib Jpeg =================================================================
#ExternalProject_Add(libjpeg
#  PREFIX ${GraphLab_SOURCE_DIR}/deps/libjpeg
#  URL http://www.ijg.org/files/jpegsrc.v8d.tar.gz
#  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
#  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
#  BUILD_IN_SOURCE 1)
#


# v8 ======================================================================

# ImageMagick =================================================================
# # set(imagemagick_configure_cmd
# #   bash -c "C_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include CPLUS_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib LDFLAGS='-L${GraphLab_SOURCE_DIR}/deps/local/lib' ./configure  --prefix=${GraphLab_SOURCE_DIR}/deps/local --without-threads --disable-openmp --disable-opencl")
# # ExternalProject_Add(imagemagick
# #   PREFIX ${GraphLab_SOURCE_DIR}/deps/imagemagick
# #   URL http://www.imagemagick.org/download/ImageMagick.tar.gz
# #   # URL_MD5 010b63a2542c4ec4918c8cb431c00356
# #   BUILD_IN_SOURCE 1
# #   CONFIGURE_COMMAND ${imagemagick_configure_cmd}
# #   BUILD_COMMAND
# #   C_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
# #   CPLUS_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
# #   LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib
# # #  LDFLAGS=-L${GraphLab_SOURCE_DIR}/deps/local/lib
# #   make
# #   INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)
#ExternalProject_Add(imagemagick
#  PREFIX ${GraphLab_SOURCE_DIR}/deps/imagemagick
#  URL http://www.imagemagick.org/download/ImageMagick.tar.gz
#  # URL_MD5 010b63a2542c4ec4918c8cb431c00356
#  CONFIGURE_COMMAND
#  C_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
#  CPLUS_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
#  LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib
#  LDFLAGS=-L${GraphLab_SOURCE_DIR}/deps/local/lib
#  <SOURCE_DIR>/configure
#  --prefix=<INSTALL_DIR>
#  --without-threads
#  --disable-openmp
#  --disable-opencl
#  BUILD_COMMAND
#  C_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
#  CPLUS_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include
#  LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib
#  LDFLAGS=-L${GraphLab_SOURCE_DIR}/deps/local/lib
#  make
#  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)
#include_directories(${GraphLab_SOURCE_DIR}/deps/local/include/ImageMagick)
#add_dependencies(imagemagick libjpeg)

# macro(requires_imagemagick NAME)
#  target_link_libraries(${NAME} Magick++ jpeg)
#  add_dependencies(${NAME} imagemagick libjpeg)
# endmacro(requires_imagemagick)

# OpenCV =================================================================
ExternalProject_Add(opencv
  PREFIX ${GraphLab_SOURCE_DIR}/deps/opencv
  URL http://sourceforge.net/projects/opencvlibrary/files/opencv-unix/2.4.9/opencv-2.4.9.zip/download
  # URL_MD5 010b63a2542c4ec4918c8cb431c00356
  PATCH_COMMAND patch -N -p0 cmake/OpenCVModule.cmake -i ${GraphLab_SOURCE_DIR}/patches/opencv_apple_rpath.patch || true
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DBUILD_TESTS=OFF -DBUILD_PERF_TESTS=OFF -DBUILD_PACKAGE=OFF -DBUILD_EXAMPLES=OFF -DOPENCV_BUILD_3RDPARTY_LIBS=ON -DBUILD_SHARED_LIBS=ON -DBUILD_DOCS=OFF -DBUILD_JPEG=ON -DCMAKE_INCLUDE_PATH=${GraphLab_SOURCE_DIR}/deps/local/include -DWITH_CUBLAS=OFF -DWITH_1394=OFF -DWITH_AVFOUNDATION=OFF -DWITH_CUDA=OFF -DWITH_CUFFT=OFF -DWITH_FFMPEG=OFF -DWITH_GSTREAMER=OFF -DWITH_GTK=OFF -DWITH_QUICKTIME=OFF -DWITH_VIDEOINPUT=OFF -DWITH_XIMEA=OFF -DWITH_XINE=OFF -DWITH_V4L=OFF -DWITH_UNICAP=OFF -DWITH_QT=OFF -DWITH_JASPER=NO -DWITH_TIFF=NO -DWITH_OPENCL=OFF -DCMAKE_LIBRARY_PATH=${GraphLab_SOURCE_DIR}/deps/local/lib
  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)

# add_dependencies(opencv eigen)

macro(requires_opencv NAME)
  get_target_property(tmp ${NAME} COMPILE_FLAGS)
  if (NOT tmp)
    set(tmp "-I${CMAKE_SOURCE_DIR}/deps/local/include/opencv")
  else()
    set(tmp "${tmp} -I${CMAKE_SOURCE_DIR}/deps/local/include/opencv")
  endif()
  set_target_properties(${NAME} PROPERTIES COMPILE_FLAGS "${tmp}")
  target_link_libraries(${NAME} opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy opencv_contrib opencv_stitching)
  add_dependencies(${NAME} opencv)
endmacro(requires_opencv)

# zlib  =======================================================================
# We rely on the zlib library to read gzip compressed files
# (using boost iostreams).
# check_library_exists(z zlibVersion "" ZLIB_FOUND)
# if (ZLIB_FOUND)
#   add_definitions(-DHAS_ZLIB)
#   link_libraries(z)
# endif()

# libevent ====================================================================
# LibEvent is used in the RPC layer to manage the interaction between the
# TCP stack and the event handler threads

if (APPLE)
  ExternalProject_Add(libevent
    PREFIX ${GraphLab_SOURCE_DIR}/deps/event
    URL http://iweb.dl.sourceforge.net/project/levent/libevent/libevent-2.0/libevent-2.0.18-stable.tar.gz
    URL_MD5 aa1ce9bc0dee7b8084f6855765f2c86a
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-openssl --enable-shared=no
    INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
  )
else()
  ExternalProject_Add(libevent
    PREFIX ${GraphLab_SOURCE_DIR}/deps/event
    URL http://iweb.dl.sourceforge.net/project/levent/libevent/libevent-2.0/libevent-2.0.18-stable.tar.gz
    URL_MD5 aa1ce9bc0dee7b8084f6855765f2c86a
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --disable-openssl --enable-shared=no
    INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
    INSTALL_COMMAND prefix=<INSTALL_DIR>/ make install && ${GraphLab_SOURCE_DIR}/patches/libevent_clean_and_remap.sh <INSTALL_DIR>/lib
  )
endif()
#link_libraries(event)
#link_libraries(event_pthreads)


# libjson ====================================================================
# Lib Json is used to support json serialization for long term storage of
# graph data.
#ExternalProject_Add(libjson
#  PREFIX ${GraphLab_SOURCE_DIR}/deps/json
#  URL http://graphlab.org/deps/libjson_7.6.0.zip
#  URL_MD5 dcb326038bd9b710b8f717580c647833
#  BUILD_IN_SOURCE 1
#  CONFIGURE_COMMAND ""
#  PATCH_COMMAND patch -N -p1 -i ${GraphLab_SOURCE_DIR}/patches/libjson.patch || true
#  BUILD_COMMAND make
#  INSTALL_COMMAND prefix=<INSTALL_DIR>/ make install
#  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
#  )


# zookeeper ================================================================

ExternalProject_Add(zookeeper
  PREFIX ${GraphLab_SOURCE_DIR}/deps/zookeeper
  URL http://mirror.metrocast.net/apache/zookeeper/zookeeper-3.5.1-alpha/zookeeper-3.5.1-alpha.tar.gz
  URL_MD5 d85f9751724d3f20f792803b61c4db24
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy_directory ${GraphLab_SOURCE_DIR}/patches/zookeeper/ <SOURCE_DIR>
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --disable-shared
  INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)

## =============================================================
# JNI =========================================================================
# We use JNI for jvm interfacing and for hdfs support.  We currently only
# detect the presence of java and do not install our own version.

if(NOT NO_JAVAC)
  find_package(Java)
  find_package(JNI)


  if((EXISTS ${JAVA_INCLUDE_PATH}/jni.h))
    set(JNI_REALLY_FOUND 1 CACHE BOOL "Tests whether the header actually exists")
    message(STATUS "jni.h was found at " ${JAVA_INCLUDE_PATH}/jni.h)
    include_directories(${JNI_INCLUDE_DIRS})
    ## Guess java home location
    if(EXISTS $ENV{JAVA_HOME} )
      set(JAVA_HOME $ENV{JAVA_HOME} CACHE STRING "Location of Java Home")
      message(STATUS "Java home set by user: " ${JAVA_HOME})
    else()
      string(REGEX REPLACE "\\/include" "" JAVA_HOME_GUESS  ${JAVA_AWT_INCLUDE_PATH})
      set(JAVA_HOME ${JAVA_HOME_GUESS} CACHE STRING "Location of Java Home")
      message(STATUS "Java Home guessed: " ${JAVA_HOME})
    endif()
  else( )
    set(JNI_REALLY_FOUND 0 CACHE BOOL "Tests whether the header actually exists")
    message(STATUS "jni.h was not found at " ${JAVA_INCLUDE_PATH}/jni.h)
  endif( )

  # Ant =========================================================================
  # Ant is used to build the graphlab JVM interface
  find_package(Ant)

  # LibHDFS =====================================================================
  # If JNI is found we install libhdfs which allows programs to read and write
  # to hdfs filesystems
  if( JNI_REALLY_FOUND )
    message(STATUS "Building libhdfs")
    ExternalProject_Add(hadoop
      PREFIX ${GraphLab_SOURCE_DIR}/deps/hadoop
      URL https://archive.apache.org/dist/hadoop/core/hadoop-1.0.1/hadoop-1.0.1.tar.gz
#      URL http://www.gtlib.gatech.edu/pub/apache/hadoop/common/hadoop-1.0.1/hadoop-1.0.1.tar.gz
      URL_MD5 e627d9b688c4de03cba8313bd0bba148
      UPDATE_COMMAND chmod +x <SOURCE_DIR>/src/c++/libhdfs/install-sh <SOURCE_DIR>/src/c++/libhdfs/configure
      PATCH_COMMAND patch -N -p1 -i ${GraphLab_SOURCE_DIR}/patches/libhdfs.patch || true
      BUILD_IN_SOURCE 1
      CONFIGURE_COMMAND <SOURCE_DIR>/src/c++/libhdfs/configure JVM_ARCH=tune=generic --prefix=<INSTALL_DIR> --with-java=${JAVA_HOME} --enable-shared=no --enable-static=yes
      INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
      TEST_COMMAND bash -c "mkdir -p <INSTALL_DIR>/include && cp <SOURCE_DIR>/src/c++/libhdfs/hdfs.h <INSTALL_DIR>/include")
    set(HADOOP_FOUND 1 CACHE BOOL "Hadoop was found")
    add_definitions(-DHAS_HADOOP)
    # link_libraries(hdfs ${JAVA_JVM_LIBRARY})
  endif( JNI_REALLY_FOUND )



endif()



# MPI =========================================================================
# MPI is currently used to establish communication between rpc nodes
#sgr this name is confusing. i like WITH_MPI better. likewise for NO_OPENMP
if(NO_MPI)
  message(STATUS "MPI disabled as a configure option.")
else()
  #find_package(MPICH2)
  find_package(MPI) 
  if(MPI_FOUND)
    message(STATUS "MPI Found: " ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    add_definitions(-DHAS_MPI)
    include_directories(${MPI_INCLUDE_PATH})
  else(MPI_FOUND)
    message(STATUS "MPI Not Found! Distributed Executables will not be compiled")
    set(MPI_LIBRARY "")
    set(MPI_EXTRA_LIBRARY "")
    # ExternalProject_Add(openmpi
    #   PREFIX ${GraphLab_SOURCE_DIR}/deps/openmpi
    #   URL http://www.open-mpi.org/software/ompi/v1.4/downloads/openmpi-1.4.5.tar.gz
    #   CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    #   INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local)
  endif(MPI_FOUND)
endif()



#ExternalProject_Add(qthreads
#                    PREFIX ${GraphLab_SOURCE_DIR}/deps/qthreads
#                    URL http://graphlab.org/deps/qthreads.tar.bz2
#                    URL_MD5 fd94ae41a06f3756c1042d47cfa671d3
#                    INSTALL_DIR ${GraphLab_SOURCE_DIR}/deps/local
#                    BUILD_IN_SOURCE 1
#                    CONFIGURE_COMMAND ./configure --prefix=<INSTALL_DIR> --enable-oversubscription --with-scheduler=sherwood --enable-lf-febs -enable-condwait-queue
#                    )



## =============================================================
# Core graphlab requirements
macro(requires_core_deps NAME)
  target_link_libraries(${NAME}
    ${Boost_LIBRARIES}
    z
    event event_pthreads
    zookeeper_mt)
  add_dependencies(${NAME} boost libevent zookeeper)
  if(MPI_FOUND)
    target_link_libraries(${NAME} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  endif(MPI_FOUND)
  if(HADOOP_FOUND)
    target_link_libraries(${NAME} hdfs ${JAVA_JVM_LIBRARY})
    add_dependencies(${NAME} hadoop)
  endif(HADOOP_FOUND)
  if(NOT NO_TCMALLOC)
    target_link_libraries(${NAME} tcmalloc)
    add_dependencies(${NAME} libtcmalloc)
  endif()
endmacro(requires_core_deps)





## ============================================================================
# System Checks

## ============================================================================
# Test for cpu affinity support
set(CMAKE_REQUIRED_LIBRARIES "pthread")
check_function_exists(pthread_setaffinity_np HAS_SET_AFFINITY)
set(CMAKE_REQUIRED_LIBRARIES ${crlbackup})

include(CheckCXXCompilerFlag)
## ============================================================================
# check if MARCH is set
if (NOT MARCH)
  #set to native if supported
  check_cxx_compiler_flag(-march=native HAS_MARCH_NATIVE)
  if(HAS_MARCH_NATIVE)
    set(MARCH "native")
  else()
    set(MARCH "x86-64")
  endif()
endif()

# check if MTUNE is set
if (NOT MTUNE)
  #set to native if supported
  check_cxx_compiler_flag(-mtune=native HAS_MTUNE_NATIVE)
  if(HAS_MTUNE_NATIVE)
    set(MTUNE "native")
  else()
    set(MTUNE "x86-64")
  endif()
endif()

## ============================================================================
# Setup compiler

# If profiling mode is enabled then additional profiling flags are set for
# the compiler
if (COMPILE_PROFILING MATCHES 1)
  set(PROFILING_FLAGS "-DUSE_EVENT_LOG -DUSE_TRACEPOINT")
else()
  set(PROFILING_FLAGS "")
endif()

# check for SSE instruction set (needed for CRC32)
#check_cxx_compiler_flag(-msse4.2 HAS_SSE42)
#set(SSE42_FLAG "")
#if (HAS_SSE42)
# set(SSE42_FLAG "-msse4.2")
#endif()


#disable Wno-unused-local-typedefs if available
check_cxx_compiler_flag(-Wno-unused-local-typedefs HAS_WNO_LOCAL_TYPEDEFS)
if(HAS_WNO_LOCAL_TYPEDEFS)
  set(WNO_LOCAL_TYPEDEFS "-Wno-unused-local-typedefs")
else()
  set(WNO_LOCAL_TYPEDEFS "")
endif()


# Set the debug flags
set(CMAKE_C_FLAGS_DEBUG
  "-O0 -Wno-attributes -march=${MARCH} -Winit-self ${PROFILING_FLAGS} ${COMPILER_FLAGS}"
  CACHE STRING "compiler options" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG
  "-O0 ${WNO_LOCAL_TYPEDEFS} -Wno-attributes -march=${MARCH} -Winit-self ${PROFILING_FLAGS} ${COMPILER_FLAGS}"
  CACHE STRING "compiler options" FORCE)

set(CMAKE_C_FLAGS_RELEASE
  "-O3 -Wno-attributes -march=${MARCH} -mtune=${MTUNE} ${PROFILING_FLAGS} ${COMPILER_FLAGS}"
  CACHE STRING "compiler options" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE
  "-O3 ${WNO_LOCAL_TYPEDEFS} -Wno-attributes -march=${MARCH} -mtune=${MTUNE} ${PROFILING_FLAGS} ${COMPILER_FLAGS}"
  CACHE STRING "compiler options" FORCE)

set(CMAKE_REQUIRED_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
check_cxx_source_compiles("int main(int argc, char** argv) { return __builtin_ia32_crc32di(0,0);}" HAS_CRC32)
set(CMAKE_REQUIRED_FLAGS "")
if(HAS_CRC32)
  set(CMAKE_CXX_FLAGS_RELEASE
      "${CMAKE_CXX_FLAGS_RELEASE} -DHAS_BUILTIN_CRC32"
      CACHE STRING "compiler options" FORCE)
endif()

# check for gcc 4.8


# # Mex setup
# set(CMAKE_CXX_FLAGS_MEX
#   "${CMAKE_CXX_FLAGS_RELEASE} -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread ")
# set(CMAKE_C_FLAGS_MEX
#   "${CMAKE_C_FLAGS_RELEASE}   -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread ")


if (CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "Release build with C++ flags: " ${CMAKE_CXX_FLAGS_RELEASE})
  message(STATUS "Release build with C flags: " ${CMAKE_C_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "Debug build with C++ flags: " ${CMAKE_CXX_FLAGS_DEBUG})
  message(STATUS "Debug build with C flags: " ${CMAKE_C_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES "Mex")
  message(STATUS "Mex CXX" ${CMAKE_CXX_FLAGS_MEX})
  message(STATUS "Mex C" ${CMAKE_C_FLAGS_MEX})
else()
  message(WARNING "Unknown build type: " ${CMAKE_BUILD_TYPE} "!")
endif()







## ============================================================================
## ============================================================================
## ============================================================================
# Setup testing tools
# Make sure testing is enabled
enable_testing()
# Use Python interpreter
find_package(PythonInterp)
set(CXXTESTGEN ${CMAKE_SOURCE_DIR}/cxxtest/cxxtestgen)
# create a macro to define a test
macro(ADD_CXXTEST NAME)
  if(PYTHONINTERP_FOUND)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp
      COMMAND
      ${PYTHON_EXECUTABLE} ${CXXTESTGEN}
      --runner=ErrorPrinter
      -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp ${ARGV}
      DEPENDS ${ARGV}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
  endif(PYTHONINTERP_FOUND)
  add_graphlab_executable(${NAME}test ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp)
  set_source_files_properties( ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp
    PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}" )

  add_test(${NAME} ${NAME}test)
endmacro(ADD_CXXTEST)


# macro(ADD_CPPTEST NAME)
#   add_graphlab_executable(${NAME}test ${NAME})
#   set_source_files_properties(${NAME}
#     PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}" )
#   add_test(${NAME} ${NAME}test)
# endmacro(ADD_CPPTEST)






## ============================================================================
## ============================================================================
## ============================================================================
# Build Macros

# copy_file ===================================================================
# copy a single file into build environment
macro(copy_file NAME)
  message(STATUS "Copying File: " ${NAME})
  file(INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/${NAME}
    DESTINATION   ${CMAKE_CURRENT_BINARY_DIR} )
endmacro(copy_file NAME)

# copy_files ==================================================================
# copy all files matching a pattern into the build environment
macro(copy_files NAME)
  message(STATUS "Copying Files: " ${NAME})
  file(INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION  ${CMAKE_CURRENT_BINARY_DIR}
    FILES_MATCHING PATTERN ${NAME} )
endmacro(copy_files NAME)


# macro(EXEC file prog)
#   add_custom_command(
#     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}
#     COMMAND
#     ${prog}
#     DEPENDS ${prog}
#     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#     )
# endmacro(EXEC)

# add_graphlab_executable =====================================================
# copy all files matching a pattern into the build environment
macro(add_graphlab_executable NAME)
  add_executable(${NAME} ${ARGN})
  target_link_libraries(${NAME} graphlab)
  add_dependencies(${NAME} graphlab)
  requires_core_deps(${NAME})
endmacro(add_graphlab_executable)


macro(add_extension_executable NAME)
  add_executable(${NAME} extension_main.cpp ${ARGN})
#  if(!APPLE)
#set_target_properties(${NAME} PROPERTIES LINK_FLAGS -Wl,-wrap,main)
#  endif()
  target_link_libraries(${NAME} graphlab graphlab_extension)
  add_dependencies(${NAME} graphlab graphlab_extension)
  requires_core_deps(${NAME})
endmacro(add_extension_executable)



add_custom_target(external_dependencies)
requires_core_deps(external_dependencies)

# add_jni_library =============================================================
# If jni was found then create a jni library.  Otherwise generate a warning
macro(add_jni_library NAME)
  # Only build if JNI was found
  if (JNI_REALLY_FOUND)
    include_directories(
      ${JAVA_INCLUDE_PATH}
      ${JAVA_INCLUDE_PATH}/linux
      )
    message(STATUS "Detected JNI library " ${NAME})
    add_library(${NAME} SHARED ${ARGN})
    target_link_libraries(${NAME}  graphlab_pic)
    target_link_libraries(${NAME} ${Boost_SHARED_LIBRARIES})
    # IF (APPLE)
    #         SET(CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS "-dynamiclib -m64 ")
    #         SET_TARGET_PROPERTIES(${NAME} PROPERTIES SUFFIX .jnilib)
    # ENDIF (APPLE)
  elseif ()
    message(STATUS "Not building " ${NAME} " because JNI was not found")
  endif ()
endmacro(add_jni_library)
