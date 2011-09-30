## Compiler Switches must come before anything else
# Force the use of gcc-4.2 as the default compiler on macs
if(APPLE)
  set(CMAKE_C_COMPILER "gcc-4.2")	
  set(CMAKE_CXX_COMPILER "c++-4.2")
endif(APPLE)




## ===================================================================
# Basic Project Declaration 

project(GraphLab)

include (CheckLibraryExists) 
include (CheckFunctionExists) 

cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0005 NEW)
endif(COMMAND cmake_policy)

# Determine where additional GraphLab specific cmake modules are
# defined
set(CMAKE_MODULE_PATH ${GraphLab_SOURCE_DIR}/cmake)




## ===================================================================
# Process command line arguments defined using cmake -D arg=val

if(NOT EXPERIMENTAL)
  set(EXPERIMENTAL 0 CACHE BOOL "Enables building of experimental
  components.")
endif()

if(NOT KC_ROOT)
  set(KC_ROOT "" CACHE STRING "Kyoto Cabinet Prefix")
endif()

# if(YRL_EXPERIMENTAL)
#   message(STATUS 
#     "\n"
#     "===============================================================\n"
#     "| Building Experimental Yahoo Labs Research Components.       |\n"
#     "|                                                             |\n"
#     "| These components are being included in the trunk to         |\n"
#     "| facilitate rapid incorporation.  The primary features       |\n"
#     "| under active development include:                           |\n"
#     "|                                                             |\n"
#     "|    a) An Hadoop interface for large graph construction.     |\n"
#     "|    b) An Hadoop interface for launching GraphLab on         |\n"
#     "|       hadoop based clouds                                   |\n"
#     "|    c) An improved sync interface that is more similar to    |\n"
#     "|       traditional map-reduce and provides more fine-grain   |\n"
#     "|       control over sync activation.                         |\n"
#     "|                                                             |\n" 
#     "===============================================================")
#   add_definitions(-DYRL_EXPERIMENTAL)
#   message(STATUS "Using Avro Path: " ${AVRO_ROOT})
#   include_directories(${AVRO_ROOT}/include)
#   link_directories(${AVRO_ROOT}/lib)
# endif()




## ===================================================================
# Setup Boost Libraries

# Use the Boost libraries
set(Boost_USE_STATIC_LIBS TRUE)
find_package(Boost 1.39 COMPONENTS 
  program_options; filesystem; system; iostreams REQUIRED)

set(Boost_SHARED_LIBRARIES "")
# Stupid hack because Cmake cannot seperate the shared from the static versions
foreach(blib ${Boost_LIBRARIES})
  message(STATUS "Boost libs: " ${blib})
  string(REGEX REPLACE "\\.a$" ${CMAKE_SHARED_LIBRARY_SUFFIX} bout ${blib})
  message(STATUS "Boost dyn libs: " ${bout})
  set(Boost_SHARED_LIBRARIES ${Boost_SHARED_LIBRARIES} ${bout})
endforeach()
message(STATUS "Boost Shared libs: " ${Boost_SHARED_LIBRARIES})

## Not used currently
# find_package(Java)


## ===================================================================
# Setup CUDA libraries
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


if (HAS_EIGEN)
  message(STATUS "Found Eigen linear algebra package")
else ()
  message(STATUS "Eigen linear algebra package not found")
endif()



#Find the JNI library
find_package(JNI)
if(EXISTS ${JAVA_INCLUDE_PATH}/jni.h  )
  set(JNI_REALLY_FOUND 1 CACHE BOOL "Tests whether the header
  actually exists")
message(STATUS "jni.h was found at " ${JAVA_INCLUDE_PATH}/jni.h)
else( )
  set(JNI_REALLY_FOUND 0 CACHE BOOL "Tests whether the header
  actually exists")
message(STATUS "jni.h was not found at " ${JAVA_INCLUDE_PATH}/jni.h)
endif( )




# EXPERIMENTAL (distributed) Capabilities
find_package(MPICH2)

if(MPI_FOUND)
  message(STATUS "MPI Found: " ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  add_definitions(-DHASMPI)
else(MPI_FOUND)
  message(STATUS "MPI Not Found! Distributed Executables will not be compiled")
endif(MPI_FOUND)


if(EXPERIMENTAL)
  #Find zoltan library
  set(ZOLTAN_ROOT "/usr/local/glzoltan" CACHE STRING 
    "The path to zoltan/include and zoltan/lib")
  if(EXISTS ${ZOLTAN_ROOT} )
    set(ZOLTAN_FOUND 1 CACHE BOOL "True if the Zoltan library was found")
    set(ZOLTAN_INCLUDE_DIR 
      ${ZOLTAN_ROOT}/include 
      CACHE STRING "Zoltan library include directory")
    set(ZOLTAN_LINK_DIR 
      ${ZOLTAN_ROOT}/lib 
      CACHE STRING "Zoltan library lib directory")
    set(ZOLTAN_LIBS 
      zoltan
      ptscotch
      ptscotcherr
      parmetis
      metis
      CACHE STRING "Zoltan library include directory")
    message(STATUS "Zoltan library was found at " ${ZOLTAN_ROOT})

    include_directories( ${ZOLTAN_INCLUDE_DIR} )
    link_directories( ${ZOLTAN_LINK_DIR}  )
  else( )
    set(ZOLTAN_FOUND 0 CACHE BOOL "True if the Zoltan library was found")
    message(STATUS "Zoltan library was NOT found at " ${ZOLTAN_ROOT})
  endif( )

  # check for SCTP
  find_library(Sctp NAMES sctp)
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



else()
  set(Sctp-FOUND)
  set(Sctp-NOTFOUND 1)
endif()



find_library(KYOTO_SHARED_LIBRARY kyotocabinet PATHS ${KC_ROOT}/lib)
find_path(KYOTO_INCLUDE_DIR kchashdb.h PATHS ${KC_ROOT}/include)
if (KYOTO_SHARED_LIBRARY MATCHES KYOTO_LIBRARY-NOTFOUND)
  message(FATAL_ERROR "Kyoto Cabinet not found. Run bootstrap!")
endif()

if (KYOTO_INCLUDE_DIR MATCHES KYOTO_INCLUDE_DIR-NOTFOUND)
  message(FATAL_ERROR "Kyoto Cabinet includes not found. Run bootstrap!")
else()
  message(STATUS "Kyoto Cabinet Found")
endif()
# get the static version
string(REGEX REPLACE "\\${CMAKE_SHARED_LIBRARY_SUFFIX}$" ".a" KYOTO_STATIC_LIBRARY ${KYOTO_SHARED_LIBRARY})
message(STATUS "Kyoto Cabinet: ${KYOTO_STATIC_LIBRARY}" )


include_directories(${KYOTO_INCLUDE_DIR})


check_library_exists(z zlibVersion "" ZLIB_FOUND)
if (ZLIB_FOUND)
  add_definitions(-DHAS_ZLIB)
  link_libraries(z)
endif()


# check for itpp
include(CheckCXXSourceCompiles)

set(crlbackup ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES "itpp")

check_cxx_source_compiles("#include <itpp/itbase.h>\nint main(int
argc, char** argv) { itpp::vec v; return 0; }" ITPP-FOUND)

set(CMAKE_REQUIRED_LIBRARIES ${crlbackup})  
if (ITPP-FOUND)
  link_libraries(itpp)
  add_definitions(-DHAS_ITPP)
  message(STATUS "ITPP found.")
else()
  message(STATUS "ITPP was not found")
endif()

if (DEFINED ITPP_INCLUDE_DIR)
  include_directories(${ITPP_INCLUDE_DIR})
  link_libraries(itpp)
  link_libraries(lapack)
  link_libraries(blas)
  add_definitions(-DHAS_ITPP)
  message(STATUS "Manually setting itpp include dir: " ${ITPP_INCLUDE_DIR})
  set(ITPP-FOUND 1)
endif() 

if (DEFINED ITPP_LINK_DIR)
  link_directories(${ITPP_LINK_DIR})
  link_libraries(lapack)
  link_libraries(blas)
  link_libraries(itpp)
  add_definitions(-DHAS_ITPP)
  message(STATUS "Manually setting itpp link dir: " ${ITPP_LINK_DIR})
  set(ITPP-FOUND 1)
endif() 



set(crlbackup ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES "pthread")
check_function_exists(pthread_setaffinity_np HAS_SET_AFFINITY) 
set(CMAKE_REQUIRED_LIBRARIES ${crlbackup})  




# test for system specific features



# check for TCMalloc
check_library_exists("tcmalloc -lpthread" malloc "" TCMALLOC-FOUND)
if (TCMALLOC-FOUND)
  message(STATUS "tcmalloc found. Using it everywhere")
  link_libraries(tcmalloc)
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
if (HAVE_64_BIT MATCHES 1)
  set(MARCH "x86-64")
else()
  set(MARCH "i686")
endif()

set(COMMON_FLAGS "-Wall -g -fopenmp")

if (APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CNAME_EXE_LINKER_FLAGS} -Xlinker -no_pie")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -O0 -march=${MARCH} -Winit-self ${COMMON_FLAGS}") 
#-fno-implicit-templates")
### Set the release flags --------------------------------------------------->
if (CROSS_COMPILE_32_BIT MATCHES 1)
  set(CMAKE_CXX_FLAGS_RELEASE 
    "-m32 -O3 -Wextra -march=i686 -mtune=i686 ${COMMON_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE 
    "-m32 -O3 -Wextra -march=i686 -mtune=i686 ${COMMON_FLAGS}")
else(CROSS_COMPILE_32_BIT MATCHES 1)
  
  # Set the C compiler flags in release mode on 64bit systems
  set(CMAKE_C_FLAGS_RELEASE 
    "-O3 -march=${MARCH} -mtune=native ${COMMON_FLAGS}")
  if (COMPILE_PROFILING MATCHES 1) 
    set(CMAKE_CXX_FLAGS_RELEASE 
      "-O2 -march=${MARCH} -mtune=native -pg ${COMMON_FLAGS}")
  else(COMPILE_PROFILING MATCHES 1)
    set(CMAKE_CXX_FLAGS_RELEASE 
      "-O3 -march=${MARCH} -mtune=native ${COMMON_FLAGS}")
  endif(COMPILE_PROFILING MATCHES 1)

endif(CROSS_COMPILE_32_BIT MATCHES 1)

# Mex setup
set(CMAKE_CXX_FLAGS_MEX  "${CMAKE_CXX_FLAGS_RELEASE} -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread ")
set(CMAKE_C_FLAGS_MEX  "${CMAKE_C_FLAGS_RELEASE}   -D_GNU_SOURCE  -fexceptions -fPIC -fno-omit-frame-pointer -pthread ")


if (CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "Release build with C++ flags: " ${CMAKE_CXX_FLAGS_RELEASE})
  message(STATUS "Release build with C flags: " ${CMAKE_C_FLAGS_RELEASE})
  #add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_RELEASE})
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "Debug build with C++ flags: " ${CMAKE_CXX_FLAGS_DEBUG})
  message(STATUS "Debug build with C flags: " ${CMAKE_C_FLAGS_DEBUG})     
  #add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_DEBUG})
elseif(CMAKE_BUILD_TYPE MATCHES "Mex")
  message(STATUS "Mex CXX" ${CMAKE_CXX_FLAGS_MEX})
  message(STATUS "Mex C" ${CMAKE_C_FLAGS_MEX})
  #add_definitions(-DCOMPILEFLAGS=${CMAKE_CXX_FLAGS_MEX})
else()
  message(WARNING "Unknown build type: " ${CMAKE_BUILD_TYPE} "!")
endif()

##endif(CMAKE_COMPILER_IS_GNUCXX)





# uncomment the following line to get CMake variables to print to screen
# include(CMakePrintSystemInformation)

## ===================================================================
# testing tools
# Make sure testing is enabled
enable_testing()
# Use Python interpreter
find_package(PythonInterp)
set(CXXTESTGEN ${CMAKE_SOURCE_DIR}/cxxtest/cxxtestgen.py)
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
  add_test(${NAME} ${NAME}test)
endmacro(ADD_CXXTEST)



## ===================================================================
# Useful macros 

macro(copyfile NAME)
  add_custom_target(
    ${NAME} ALL
    COMMAND 
    ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${NAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${NAME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Copying ${NAME}"
    )
endmacro(copyfile NAME)

macro(EXEC file prog)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}
    COMMAND
    ${prog}
    DEPENDS ${prog}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
endmacro(EXEC)


macro(add_graphlab_executable NAME)
  add_executable(${NAME} ${ARGN})
  target_link_libraries(${NAME} ${Boost_LIBRARIES})
  target_link_libraries(${NAME} ${GraphLab_LIBRARIES})
  target_link_libraries(${NAME} ${KYOTO_STATIC_LIBRARY})
endmacro(add_graphlab_executable)





macro(add_dist2_executable NAME)
  if(MPI_FOUND)
    add_executable(${NAME} ${ARGN})
    target_link_libraries(${NAME} 
      graphlab_dist2 
      graphlab
      ${MPI_LIBRARY} 
      ${MPI_EXTRA_LIBRARY} 
      ${Boost_LIBRARIES}
      ${KYOTO_STATIC_LIBRARY}
      )
    set_source_files_properties( ${ARGN}
      PROPERTIES COMPILE_DEFINITIONS "GLDISTRIBUTED2" )
    # set_property(TARGET ${NAME} APPEND PROPERTY 
    #   COMPILE_FLAGS -I${MPI_INCLUDE_PATH})
    #    add_definitions(-DGLDISTRIBUTED)

  else(MPI_FOUND)
    message(STATUS "Not Building " ${NAME} " since mpi was not found.")
  endif(MPI_FOUND)
endmacro(add_dist2_executable)






macro(add_zoltan_executable NAME)
  if(MPI_FOUND AND ZOLTAN_FOUND)
    add_executable(${NAME} ${ARGN})   
    target_link_libraries(${NAME} 
      graphlab
      graphlab_dist2
      graphlab_zoltan
      ${ZOLTAN_LIBS}
      ${MPI_LIBRARY} 
      ${MPI_EXTRA_LIBRARY} )
  else( )
    if(EXPERIMENTAL)   
      message(STATUS "Not Building " ${NAME} " because zoltan was not found")
    endif()
  endif( )
endmacro(add_zoltan_executable)



MACRO(dist2_target_link_libraries NAME)
  IF(MPI_FOUND)
    TARGET_LINK_LIBRARIES(${NAME} ${ARGN})
  ENDIF(MPI_FOUND)
ENDMACRO(dist2_target_link_libraries)

MACRO(dist2_set_target_properties NAME)
  IF(MPI_FOUND)
    set_target_properties(${NAME} ${ARGN})
  ENDIF(MPI_FOUND)
ENDMACRO(dist2_set_target_properties)



MACRO(add_jni_library NAME)
  # Only build if JNI was found
  IF (JNI_REALLY_FOUND)
    include_directories(
      ${JAVA_INCLUDE_PATH}
      ${JAVA_INCLUDE_PATH}/linux
      )
    ADD_LIBRARY(${NAME} MODULE ${ARGN})
    TARGET_LINK_LIBRARIES(${NAME}  graphlab_pic) 
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


## ===================================================================
# Installation rules
install(DIRECTORY src/  
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp"  PATTERN "*.h" 
  PATTERN ".svn" EXCLUDE
  )

install(CODE "execute_process(COMMAND ./scripts/install_deps.sh
  ${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR})")



