if (CMAKE_BUILD_TYPE MATCHES "Release")
  set(METIS_BUILD_OPTS "${CMAKE_C_FLAGS_RELEASE} -w -DLINUX -std=c99 -D_FILE_OFFSET_BITS=64 -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis/GKlib -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis/libmetis -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis") 
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(METIS_BUILD_OPTS "${CMAKE_C_FLAGS_DEBUG} -w -DLINUX -std=c99 -D_FILE_OFFSET_BITS=64 -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis/GKlib -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis/libmetis -I${CMAKE_CURRENT_SOURCE_DIR}/extern/metis")
endif()

IF(APPLE)
  set(METIS_BUILD_OPTS "${METIS_BUILD_OPTS} -DDARWIN")
ENDIF(APPLE)

set (metis_sources
extern/metis/GKlib/b64.c
extern/metis/GKlib/dlmalloc.c
extern/metis/GKlib/getopt.c
extern/metis/GKlib/memory.c
extern/metis/GKlib/pqueue.c
extern/metis/GKlib/string.c
extern/metis/GKlib/util.c
extern/metis/GKlib/blas.c
extern/metis/GKlib/error.c
extern/metis/GKlib/htable.c
extern/metis/GKlib/seq.c
extern/metis/GKlib/timers.c
extern/metis/GKlib/dfkvkselect.c
extern/metis/GKlib/fs.c
extern/metis/GKlib/io.c
extern/metis/GKlib/pdb.c
extern/metis/GKlib/sort.c
extern/metis/GKlib/tokenizer.c
extern/metis/libmetis/balance.c
extern/metis/libmetis/graph.c
extern/metis/libmetis/mcoarsen.c
extern/metis/libmetis/mmatch.c
extern/metis/libmetis/refine.c
extern/metis/libmetis/bucketsort.c
extern/metis/libmetis/initpart.c
extern/metis/libmetis/memory.c
extern/metis/libmetis/mmd.c
extern/metis/libmetis/rkmetis.c
extern/metis/libmetis/ccgraph.c
extern/metis/libmetis/kfmetis.c
extern/metis/libmetis/mesh.c
extern/metis/libmetis/mpmetis.c
extern/metis/libmetis/separator.c
extern/metis/libmetis/checkgraph.c
extern/metis/libmetis/kmetis.c
extern/metis/libmetis/meshpart.c
extern/metis/libmetis/mrefine2.c
extern/metis/libmetis/sfm.c
extern/metis/libmetis/cmetis.c
extern/metis/libmetis/kvmetis.c
extern/metis/libmetis/mfm2.c
extern/metis/libmetis/mrefine.c
extern/metis/libmetis/srefine.c
extern/metis/libmetis/coarsen.c
extern/metis/libmetis/kwayfm.c
extern/metis/libmetis/mfm.c
extern/metis/libmetis/mrkmetis.c
extern/metis/libmetis/stat.c
extern/metis/libmetis/compress.c
extern/metis/libmetis/kwayrefine.c
extern/metis/libmetis/mincover.c
extern/metis/libmetis/mutil.c
extern/metis/libmetis/streamio.c
extern/metis/libmetis/debug.c
extern/metis/libmetis/kwayvolfm.c
extern/metis/libmetis/minitpart2.c
extern/metis/libmetis/myqsort.c
extern/metis/libmetis/subdomains.c
extern/metis/libmetis/estmem.c
extern/metis/libmetis/kwayvolrefine.c
extern/metis/libmetis/minitpart.c
extern/metis/libmetis/ometis.c
extern/metis/libmetis/timing.c
extern/metis/libmetis/fm.c
extern/metis/libmetis/match.c
extern/metis/libmetis/mkmetis.c
extern/metis/libmetis/parmetis.c
extern/metis/libmetis/util.c
extern/metis/libmetis/fortran.c
extern/metis/libmetis/mbalance2.c
extern/metis/libmetis/mkwayfmh.c
extern/metis/libmetis/pmetis.c
extern/metis/libmetis/frename.c
extern/metis/libmetis/mbalance.c
extern/metis/libmetis/mkwayrefine.c
extern/metis/libmetis/pqueue.c)

set_property(SOURCE extern/metis/GKlib/b64.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/dlmalloc.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/getopt.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/memory.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/pqueue.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/string.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/util.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/blas.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/error.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/htable.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/seq.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/timers.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/dfkvkselect.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/fs.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/io.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/pdb.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/sort.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/GKlib/tokenizer.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/balance.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/graph.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mcoarsen.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mmatch.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/refine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/bucketsort.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/initpart.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/memory.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mmd.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/rkmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/ccgraph.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kfmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mesh.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mpmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/separator.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/checkgraph.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/meshpart.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mrefine2.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/sfm.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/cmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kvmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mfm2.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mrefine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/srefine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/coarsen.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kwayfm.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mfm.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mrkmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/stat.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/compress.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kwayrefine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mincover.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mutil.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/streamio.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/debug.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kwayvolfm.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/minitpart2.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/myqsort.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/subdomains.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/estmem.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/kwayvolrefine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/minitpart.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/ometis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/timing.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/fm.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/match.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mkmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/parmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/util.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/fortran.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mbalance2.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mkwayfmh.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/pmetis.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/frename.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mbalance.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/mkwayrefine.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
set_property(SOURCE extern/metis/libmetis/pqueue.c PROPERTY COMPILE_FLAGS "${METIS_BUILD_OPTS}")
