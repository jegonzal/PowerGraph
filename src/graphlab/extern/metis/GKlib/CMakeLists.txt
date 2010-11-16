project(GraphLab)



add_library(graphlab_GKlib STATIC
  b64.c          dlmalloc.c  getopt.c  memory.c  pqueue.c  string.c     util.c
  blas.c         error.c     htable.c  omp.c     seq.c     timers.c
  dfkvkselect.c  fs.c        io.c      pdb.c     sort.c    tokenizer.c
)

set_target_properties(graphlab_GKlib PROPERTIES COMPILE_FLAGS "${METIS_BUILD_OPTS}")


get_target_property(graphlab_GKlib_sources graphlab_GKlib SOURCES)
add_library(graphlab_GKlib_pic STATIC ${graphlab_GKlib_sources})
set_target_properties(graphlab_GKlib_pic PROPERTIES COMPILE_FLAGS "${METIS_BUILD_OPTS} -fPIC")


INSTALL(TARGETS graphlab_GKlib DESTINATION lib)



