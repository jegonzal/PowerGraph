
project(GraphLab)

add_library(graphlab_metis STATIC
  balance.c     graph.c          mcoarsen.c     mmatch.c    refine.c
  bucketsort.c  initpart.c       memory.c       mmd.c       rkmetis.c
  ccgraph.c     kfmetis.c        mesh.c         mpmetis.c   separator.c
  checkgraph.c  kmetis.c         meshpart.c     mrefine2.c  sfm.c
  cmetis.c      kvmetis.c        mfm2.c         mrefine.c   srefine.c
  coarsen.c     kwayfm.c         mfm.c          mrkmetis.c  stat.c
  compress.c    kwayrefine.c     mincover.c     mutil.c     streamio.c
  debug.c       kwayvolfm.c      minitpart2.c   myqsort.c   subdomains.c
  estmem.c      kwayvolrefine.c  minitpart.c    ometis.c    timing.c
  fm.c          match.c          mkmetis.c      parmetis.c  util.c
  fortran.c     mbalance2.c      mkwayfmh.c     pmetis.c
  frename.c     mbalance.c       mkwayrefine.c  pqueue.c
)

# add_library(metis2_shared SHARED
#   balance.c     graph.c          mcoarsen.c     mmatch.c    refine.c
#   bucketsort.c  initpart.c       memory.c       mmd.c       rkmetis.c
#   ccgraph.c     kfmetis.c        mesh.c         mpmetis.c   separator.c
#   checkgraph.c  kmetis.c         meshpart.c     mrefine2.c  sfm.c
#   cmetis.c      kvmetis.c        mfm2.c         mrefine.c   srefine.c
#   coarsen.c     kwayfm.c         mfm.c          mrkmetis.c  stat.c
#   compress.c    kwayrefine.c     mincover.c     mutil.c     streamio.c
#   debug.c       kwayvolfm.c      minitpart2.c   myqsort.c   subdomains.c
#   estmem.c      kwayvolrefine.c  minitpart.c    ometis.c    timing.c
#   fm.c          match.c          mkmetis.c      parmetis.c  util.c
#   fortran.c     mbalance2.c      mkwayfmh.c     pmetis.c
#   frename.c     mbalance.c       mkwayrefine.c  pqueue.c
# )




# SET_TARGET_PROPERTIES(metis2_shared 
#   PROPERTIES
#   OUTPUT_NAME metis2
# )

set_target_properties(graphlab_metis PROPERTIES COMPILE_FLAGS "${METIS_BUILD_OPTS}")
target_link_libraries(graphlab_metis graphlab_GKlib)

get_target_property(graphlab_metis_sources graphlab_metis SOURCES)
add_library(graphlab_metis_pic STATIC ${graphlab_metis_sources})
set_target_properties(graphlab_metis_pic PROPERTIES COMPILE_FLAGS "${METIS_BUILD_OPTS} -fPIC")
target_link_libraries(graphlab_metis_pic graphlab_GKlib_pic)

INSTALL(TARGETS graphlab_metis graphlab_metis_pic DESTINATION lib)
