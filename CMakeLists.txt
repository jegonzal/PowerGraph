# Force the use of gcc-4.2 as the default compiler on macs
IF(APPLE)
  SET(CMAKE_C_COMPILER "gcc-4.2")	
  SET(CMAKE_CXX_COMPILER "c++-4.2")
ENDIF(APPLE)


project(GraphLab)

INCLUDE (CheckLibraryExists) 
cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

# cmake_policy(SET CMP0014 OLD)




set(CMAKE_MODULE_PATH ${GraphLab_SOURCE_DIR}/cmake)

#set(CROSS_COMPILE_32_BIT 1)

# Make sure testing is enabled
ENABLE_TESTING()


# Use Python interpreter
FIND_PACKAGE(PythonInterp)

# Use the Boost libraries
SET (Boost_USE_STATIC_LIBS TRUE)
FIND_PACKAGE(Boost 1.39 COMPONENTS 
  program_options; regex; filesystem; system REQUIRED)



SET(Boost_SHARED_LIBRARIES "")
# Stupid hack because Cmake cannot seperate the shared from the static versions
foreach(blib ${Boost_LIBRARIES})
  message(STATUS "Boost libs: " ${blib})
  STRING(REGEX REPLACE "\\.a$" ${CMAKE_SHARED_LIBRARY_SUFFIX} bout ${blib})
  message(STATUS "Boost dyn libs: " ${bout})
  SET(Boost_SHARED_LIBRARIES ${Boost_SHARED_LIBRARIES} ${bout})
endforeach()
message(STATUS "Boost Shared libs: " ${Boost_SHARED_LIBRARIES})

## Not used currently
# find_package(Java)



# # Use CUDA libraries (required only for GPU GraphLab)
# find_package(CUDA)
# # If CUDA libraries were not found, ignore commands:
# #   cuda_add_executable
# IF(CUDA_FOUND)
#   message(STATUS "CUDA Found: " ${CUDA_VERSION})
#   link_libraries(${CUDA_CUBLAS_LIBRARIES})
# ELSE(CUDA_FOUND)
#   message(STATUS "CUDA Not Found! GPU Executables will not be compiled.")
#   macro(CUDA_ADD_EXECUTABLE cuda_target)
#   endmacro(CUDA_ADD_EXECUTABLE cuda_target)
# ENDIF(CUDA_FOUND)



find_package(MPICH2)
IF(MPI_FOUND)
  message(STATUS "MPI Found: " ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
ELSE(MPI_FOUND)
  message(STATUS "MPI Not Found! Distributed Executables will not be compiled")
ENDIF(MPI_FOUND)





#Find the JNI library
FIND_PACKAGE(JNI)
IF(EXISTS ${JAVA_INCLUDE_PATH}/jni.h  )
  SET(JNI_REALLY_FOUND 1 CACHE BOOL "Tests whether the header
  actually exists")
message(STATUS "jni.h was found at " ${JAVA_INCLUDE_PATH}/jni.h)
ELSE( )
  SET(JNI_REALLY_FOUND 0 CACHE BOOL "Tests whether the header
  actually exists")
message(STATUS "jni.h was not found at " ${JAVA_INCLUDE_PATH}/jni.h)
ENDIF( )






#Find zoltan library
SET(ZOLTAN_ROOT "/usr/local/glzoltan" CACHE STRING 
  "The path to zoltan/include and zoltan/lib")
IF(EXISTS ${ZOLTAN_ROOT} )
  SET(ZOLTAN_FOUND 1 CACHE BOOL "True if the Zoltan library was found")
  SET(ZOLTAN_INCLUDE_DIR 
    ${ZOLTAN_ROOT}/include 
    CACHE STRING "Zoltan library include directory")
  SET(ZOLTAN_LINK_DIR 
    ${ZOLTAN_ROOT}/lib 
    CACHE STRING "Zoltan library lib directory")
  SET(ZOLTAN_LIBS 
    zoltan
    ptscotch
    ptscotcherr
    parmetis
    metis
    CACHE STRING "Zoltan library include directory")
  message(STATUS "Zoltan library was found at " ${ZOLTAN_ROOT})

  include_directories( ${ZOLTAN_INCLUDE_DIR} )
  link_directories( ${ZOLTAN_LINK_DIR}  )
ELSE( )
  SET(ZOLTAN_FOUND 0 CACHE BOOL "True if the Zoltan library was found")
  message(STATUS "Zoltan library was NOT found at " ${ZOLTAN_ROOT})
ENDIF( )



# check for SCTP
FIND_LIBRARY(Sctp NAMES sctp)
if(Sctp MATCHES Sctp-NOTFOUND)
  message(STATUS "libsctp not found. RPC over SCTP support will not be compiled")
  set(Sctp-FOUND)
  set(Sctp-NOTFOUND 1)
else()
  add_definitions(-DHAS_SCTP)
  set(Sctp-FOUND 1)
  set(Sctp-NOTFOUND)
  message(STATUS "libsctp found")
endif()


CHECK_LIBRARY_EXISTS(z zlibVersion "" ZLIB_FOUND)
if (ZLIB_FOUND)
  add_definitions(-DZLIB_FOUND)
  link_libraries(z)
endif()


# check for itpp
include(CheckCXXSourceCompiles)

set(crlbackup ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES "itpp")
CHECK_CXX_SOURCE_COMPILES("#include <itpp/itbase.h>\nint main(int argc ,char** argv) { itpp::vec v; return 0; }" ITPP-FOUND)
set(CMAKE_REQUIRED_LIBRARIES ${crlbackup})  
if (ITPP-FOUND)
  link_libraries(itpp)
  message(STATUS "ITPP found. Now we can make beautiful matrices together.")
else()
  message(STATUS "Bummer, ITPP was not found but we don't need it anyways.")
endif()




# check for TCMalloc
CHECK_LIBRARY_EXISTS("tcmalloc -lpthread" malloc "" TCMALLOC-FOUND)
if (TCMALLOC-FOUND)
  message(STATUS "tcmalloc found. Using it everywhere")
#  link_libraries(tcmalloc)
else()
  message(STATUS "tcmalloc was not found but is optional.")
endif()


# architecture detection
# check 64 bit
if(CMAKE_SIZEOF_VOID_P MATCHES 4)
  set( HAVE_64_BIT 0 )
else(CMAKE_SIZEOF_VOID_P MATCHES 4)
  set( HAVE_64_BIT 1 )
endif(CMAKE_SIZEOF_VOID_P MATCHES 4)

if(CROSS_COMPILE_32_BIT MATCHES 1)
  set(HAVE_64_BIT 0)
endif(CROSS_COMPILE_32_BIT MATCHES 1)


# add random environment variable definitions
# deal with Boost time shtuff: we want nanoseconds!
# add_definitions(-DBOOST_DATE_TIME_POSIX_TIME_STD_CONFIG)
# add_definitions(-DBOOST_ALL_DYN_LINK)

# add_definitions(-D_SCL_SECURE_NO_WARNINGS) # disable checked_iterator warnings
# add_definitions(-D_CRT_SECURE_NO_WARNINGS) # disable fopen warnings
# add_definitions(-D_SECURE_SCL=0)  # necessary for release builds of MSVC

# Remove tracing by default
# add_definitions(-DTRACE_OFF)

## Here we use our cmake tools
#set(CMAKE_MODULE_PATH ${PRL_SOURCE_DIR}/cmake)

set(BUILD_SHARED_LIBS ON)

# set include path for this and all subdirs
include_directories(
  ${GraphLab_SOURCE_DIR}/src
  ${GraphLab_SOURCE_DIR}/cxxtest
  ${Boost_INCLUDE_DIR}
  )

IF(MPI_FOUND)
  include_directories(
    ${MPI_INCLUDE_PATH}
    )
ENDIF(MPI_FOUND)






# add boost to link directories and link libraries
link_directories(
  ${Boost_LIBRARY_DIRS}
  )


# Add posix threads on non-windows systems
if(NOT ${CMAKE_SYSTEM} MATCHES Windows)
  link_libraries(pthread gomp)   
endif(NOT ${CMAKE_SYSTEM} MATCHES Windows)




# Flags for the debug and release mode
## IMPORTANT: we currently assume a GCC compiler when setting flags
# if(CMAKE_COMPILER_IS_GNUCXX)

### Set the debug flags ----------------------------------------------------->
# Note: -fkeep-inline-functions significantly slows down the compilation
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0 -Wall -Winit-self -fopenmp") 
#-fno-implicit-templates")
### Set the release flags --------------------------------------------------->
if (CROSS_COMPILE_32_BIT MATCHES 1)
  set(CMAKE_CXX_FLAGS_RELEASE 
    "-m32 -O3 -Wall -Wextra -g -mfpmath=sse -msse2 -funroll-loops -fprefetch-loop-arrays ")
  set(CMAKE_C_FLAGS_RELEASE 
    "-m32 -O3 -Wall -Wextra -g -mfpmath=sse -msse2 -funroll-loops -fprefetch-loop-arrays ")
else(CROSS_COMPILE_32_BIT MATCHES 1)
  
  # Set the C compiler flags in release mode on 64bit systems
  set(CMAKE_C_FLAGS_RELEASE 
    "-O3 -Wall -g -mfpmath=sse -msse2 -funroll-loops -fprefetch-loop-arrays -fopenmp")
  if (COMPILE_PROFILING MATCHES 1) 
    set(CMAKE_CXX_FLAGS_RELEASE 
      "-O2 -Wall -g -mfpmath=sse -msse2 -fprefetch-loop-arrays -pg -fopenmp")
  else(COMPILE_PROFILING MATCHES 1)
    set(CMAKE_CXX_FLAGS_RELEASE 
      "-O3 -Wall -g -mfpmath=sse -msse2 -funroll-loops -fprefetch-loop-arrays -fopenmp")
  endif(COMPILE_PROFILING MATCHES 1)

endif(CROSS_COMPILE_32_BIT MATCHES 1)

# Mex setup
set(CMAKE_CXX_FLAGS_MEX  "${CMAKE_CXX_FLAGS_RELEASE} -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread ")
set(CMAKE_C_FLAGS_MEX  "${CMAKE_C_FLAGS_RELEASE}   -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread ")


if (CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "Release build with C++ flags: " ${CMAKE_CXX_FLAGS_RELEASE})
  message(STATUS "Release build with C flags: " ${CMAKE_C_FLAGS_RELEASE})
  add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "Debug build with C++ flags: " ${CMAKE_CXX_FLAGS_DEBUG})
  message(STATUS "Debug build with C flags: " ${CMAKE_C_FLAGS_DEBUG})     
  add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES "Mex")
  message(STATUS "Mex CXX" ${CMAKE_CXX_FLAGS_MEX})
  message(STATUS "Mex C" ${CMAKE_C_FLAGS_MEX})
  add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_MEX})
else()
  message(WARNING "Unknown build type: " ${CMAKE_BUILD_TYPE} "!")
endif()

##endif(CMAKE_COMPILER_IS_GNUCXX)





# uncomment the following line to get CMake variables to print to screen
# include(CMakePrintSystemInformation)


SET(CXXTESTGEN ${CMAKE_SOURCE_DIR}/cxxtest/cxxtestgen.py)

MACRO(ADD_CXXTEST NAME)
  IF(PYTHONINTERP_FOUND)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp
      COMMAND
      ${PYTHON_EXECUTABLE} ${CXXTESTGEN}
      --runner=ErrorPrinter
      -o ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp ${ARGV}
      DEPENDS ${ARGV}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
  ENDIF(PYTHONINTERP_FOUND)

  ADD_EXECUTABLE(${NAME}test ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.cpp)

  ADD_TEST(${NAME} ${NAME}test)
ENDMACRO(ADD_CXXTEST)


MACRO(EXEC file prog)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}
    COMMAND
    ${prog}
    DEPENDS ${prog}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
ENDMACRO(EXEC)


MACRO(add_dist_executable NAME)
  IF(MPI_FOUND)
    ADD_EXECUTABLE(${NAME} ${ARGN})
    TARGET_LINK_LIBRARIES(${NAME} 
      graphlab_dist 
      ${MPI_LIBRARY} 
      ${MPI_EXTRA_LIBRARY} 
      graphlab)

    SET_SOURCE_FILES_PROPERTIES( ${ARGN}
      PROPERTIES COMPILE_DEFINITIONS "GLDISTRIBUTED" )
    # set_property(TARGET ${NAME} APPEND PROPERTY 
    #   COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
    #    add_definitions(-DGLDISTRIBUTED)

  ELSE(MPI_FOUND)
    message(STATUS "Not Building " ${NAME} " since mpi was not found.")
  ENDIF(MPI_FOUND)
ENDMACRO(add_dist_executable)


MACRO(add_dist2_executable NAME)
  IF(MPI_FOUND)
    ADD_EXECUTABLE(${NAME} ${ARGN})
    TARGET_LINK_LIBRARIES(${NAME} 
      graphlab_dist2 
      ${MPI_LIBRARY} 
      ${MPI_EXTRA_LIBRARY} 
      graphlab)

    SET_SOURCE_FILES_PROPERTIES( ${ARGN}
      PROPERTIES COMPILE_DEFINITIONS "GLDISTRIBUTED2" )
    # set_property(TARGET ${NAME} APPEND PROPERTY 
    #   COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
    #    add_definitions(-DGLDISTRIBUTED)

  ELSE(MPI_FOUND)
    message(STATUS "Not Building " ${NAME} " since mpi was not found.")
  ENDIF(MPI_FOUND)
ENDMACRO(add_dist2_executable)


MACRO(add_zoltan_executable NAME)
  IF(MPI_FOUND AND ZOLTAN_FOUND)
    ADD_EXECUTABLE(${NAME} ${ARGN})   
    TARGET_LINK_LIBRARIES(${NAME} 
      graphlab
      graphlab_dist2
      graphlab_zoltan
      ${ZOLTAN_LIBS}
      ${MPI_LIBRARY} 
      ${MPI_EXTRA_LIBRARY} )
  ELSE( )
    message(STATUS "Not Building " ${NAME} 
      " because zoltan was not found")
  ENDIF( )
ENDMACRO(add_zoltan_executable)






MACRO(dist_target_link_libraries NAME)
  IF(MPI_FOUND)
    TARGET_LINK_LIBRARIES(${NAME} ${ARGN})
  ENDIF(MPI_FOUND)
ENDMACRO(dist_target_link_libraries)



MACRO(add_jni_library NAME)
  # Only build if JNI was found
  IF (JNI_REALLY_FOUND)
    include_directories(
      ${JAVA_INCLUDE_PATH}
      ${JAVA_INCLUDE_PATH}/linux
      )
    ADD_LIBRARY(${NAME} MODULE ${ARGN})
    TARGET_LINK_LIBRARIES(${NAME}  graphlab_pic) 
    TARGET_LINK_LIBRARIES(${NAME}  graphlab_util_pic) 
    TARGET_LINK_LIBRARIES(${NAME}  graphlab_metis_pic) 
    TARGET_LINK_LIBRARIES(${NAME}  graphlab_GKlib_pic) 
    TARGET_LINK_LIBRARIES(${NAME} ${Boost_SHARED_LIBRARIES})
    # IF (APPLE)
    #         SET(CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS "-dynamiclib -m64 ")
    #         SET_TARGET_PROPERTIES(${NAME} PROPERTIES SUFFIX .jnilib)
    # ENDIF (APPLE)
  ENDIF ()
ENDMACRO(add_jni_library)


# FIND_PACKAGE(Matlab)
# MACRO(add_mex_library NAME)
# IF(MATLAB_FOUND)
#     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -DMX_COMPAT_32 -DMATLAB_MEX_FILE")
#     SET(MEX_SAVE_PATH "${CMAKE_SOURCE_DIR}")

#     INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
#     ADD_LIBRARY(${NAME} STATIC ${ARGN}) 
#     TARGET_LINK_LIBRARIES(${NAME}  graphlab) 

#      ADD_CUSTOM_COMMAND(TARGET ${NAME}
#       POST_BUILD
#       COMMAND ${CMAKE_MODULE_PATH}/mex_link.sh ${NAME} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} "lib${NAME}.a ${Boost_LIBRARIES} ${CMAKE_SOURCE_DIR}/matlabmex/src/graphlab/libgraphlab.a  ${CMAKE_SOURCE_DIR}/matlabmex/extern/judy/libJudy.a ${CMAKE_SOURCE_DIR}/matlabmex/extern/include/metis2/libmetis/libmetis2.a   ${CMAKE_SOURCE_DIR}/matlabmex/extern/include/metis2/GKlib/libGKlib.a  "
#     )
#    # set(CMAKE_C_FLAGS_RELEASE  "${CMAKE_C_FLAGS_RELEASE} -fPIC")  
#    # set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -fPIC")  

#   # INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIR})
#   # ADD_LIBRARY(${NAME} MODULE ${ARGN})
#   # TARGET_LINK_LIBRARIES(${NAME} ${MATLAB_LIBRARIES} graphlab) 
#   # SET(MEX_SAVE_PATH "${CMAKE_SOURCE_DIR}")
#   # SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC -std=c99 -DMX_COMPAT_32 -DMATLAB_MEX_FILE")
#   # SET(CMAKE_CPP_FLAGS "${CMAKE_CPP_FLAGS} -fPIC -std=c99 -DMX_COMPAT_32 -DMATLAB_MEX_FILE")

# ENDIF(MATLAB_FOUND)
# ENDMACRO(add_mex_library)







set(GraphLab_LIBRARIES graphlab)


# Specify the subdirectories for this root directory This will be
# used for the class branch, not the lab branch: subdirs(src examples
# EXCLUDE_FROM_ALL tests timings projects) (Remove EXCLUDE_FROM_ALL
# for lab branch.)
subdirs(src tests)


if (IS_DIRECTORY ${GraphLab_SOURCE_DIR}/apps)
  subdirs(apps)
endif()
if (IS_DIRECTORY ${GraphLab_SOURCE_DIR}/demoapps)
  subdirs(demoapps)
endif()

subdirs(tools)

## Installation rules
## ============================================================
INSTALL(DIRECTORY src/  
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp"  PATTERN "*.h" 
  PATTERN ".svn" EXCLUDE
  )


# Install target graphlab moved to after definition


