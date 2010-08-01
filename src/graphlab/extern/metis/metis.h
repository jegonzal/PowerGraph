/*
 * metis.h 
 *
 * This file contains function prototypes and constant definitions 
 * for METIS
 *
 * Started 8/9/02
 * George
 *
 */

#ifndef METIS_H
#define METIS_H 1


/****************************************************************************
* A set of defines that can be modified by the user
*****************************************************************************/
/*--------------------------------------------------------------------------
 Specifies the width of the elementary data type that will hold information
 about vertices and their adjacency lists.

 Possible values:
   32 : Use 32 bit signed integers
   64 : Use 64 bit signed integers

 A width of 64 should be specified if the number of vertices or the total
 number of edges in the graph exceed the limits of a 32 bit signed integer
 i.e., 2^31-1.
 Proper use of 64 bit integers requires that the c99 standard datatypes
 int32_t and int64_t are supported by the compiler.
 GCC does provides these definitions in stdint.h, but it may require some
 modifications on other architectures.
--------------------------------------------------------------------------*/
#define IDXTYPEWIDTH 64


/*--------------------------------------------------------------------------
 Specifies if the __thread storage directive is available by the compiler
 to indicate thread local storage. This storage directive is available in
 most systems using gcc compiler but it may not be available in other
 systems.

 Possible values:
  0 : Not available and do not use thread local storage
  1 : It is available and the __thread modifier will be used
--------------------------------------------------------------------------*/
#define HAVE_THREADLOCALSTORAGE 0


/****************************************************************************
* In principle, nothing needs to be changed beyond this point, unless the
* int32_t and int64_t cannot be found in the normal places.
*****************************************************************************/


/* Uniform definitions for various compilers */
#if defined(_MSC_VER)
  #define COMPILER_MSC
#endif
#if defined(__ICC)
  #define COMPILER_ICC
#endif
#if defined(__GNUC__)
  #define COMPILER_GCC
#endif



#if defined(COMPILER_MSC)
  #include <ctrdefs.h>
  #define __thread __declspec( thread )

  typedef __int32                 int32_t;
  typedef __int64                 int64_t;
  typedef unsigned __int32        uint32_t;
  typedef unsigned __int64        uint64_t;
#else
  #include <stdint.h>
  #include <inttypes.h>
  #include <sys/types.h>
#endif



/*------------------------------------------------------------------------
* Undefine the following #define in order to use short idxtype as the idxtype 
*-------------------------------------------------------------------------*/
#if IDXTYPEWIDTH == 32
  #define SCNIDX  SCNd32
  #define PRIIDX  PRId32

  typedef int32_t idxtype;
#elif IDXTYPEWIDTH == 64
  #define SCNIDX  SCNd64
  #define PRIIDX  PRId64

  typedef int64_t idxtype;
#else
  #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif





/*------------------------------------------------------------------------
* Function prototypes 
*-------------------------------------------------------------------------*/

#if !defined(__cdecl)
  #define __cdecl
#endif



#ifdef __cplusplus
extern "C" {
#endif

void __cdecl METIS_EstimateMemory(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, 
                   idxtype *optype, idxtype *nbytes);

void __cdecl METIS_PartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part); 

void __cdecl METIS_WPartGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *edgecut, idxtype *part); 

void __cdecl METIS_PartGraphVKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *volume, idxtype *part);

void __cdecl METIS_WPartGraphVKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *volume, idxtype *part);

idxtype  __cdecl METIS_MeshToDualCount(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, idxtype *etype,
                   idxtype *numflag);

void __cdecl METIS_MeshToDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, idxtype *etype, 
                   idxtype *numflag,  idxtype *dxadj, idxtype *dadjncy);

idxtype  __cdecl METIS_MixedMeshToDualCount(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype * elms, 
                   idxtype *etype, idxtype *numflag, idxtype *conmat, idxtype custom);

void __cdecl METIS_MixedMeshToDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *elms, 
                   idxtype *etype, idxtype *numflag,idxtype *dxadj, idxtype *dadjncy,idxtype *conmat,
                   idxtype custom);

void __cdecl METIS_MeshToNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, 
                   idxtype *dxadj, idxtype *dadjncy);

void __cdecl METIS_MixedMeshToNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype,
                   idxtype *numflag, idxtype *dxadj, idxtype *dadjncy);

void __cdecl METIS_PartMeshNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag,
                   idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart);

void __cdecl METIS_PartMixedMeshNodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag,
                   idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart);

void __cdecl METIS_PartMeshDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, 
                   idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart, idxtype wgtflag, 
                   idxtype * vwgt);

void __cdecl METIS_PartMixedMeshDual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag,  
                   idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart, idxtype *conmat, 
                   idxtype custom, idxtype wgtflag, idxtype *vwgt);

void __cdecl METIS_mCPartGraphKway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, 
                   idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, 
                   float *rubvec, idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_mCPartGraphRecursive(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy,
                   idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts,
                   idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_mCHPartGraphRecursive(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy,
                   idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts,
                   float *ubvec, idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_mCPartGraphRecursiveInternal(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, 
                   idxtype *adjncy, float *nvwgt, idxtype *adjwgt, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part);

void __cdecl METIS_mCHPartGraphRecursiveInternal(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, 
                   idxtype *adjncy, float *nvwgt, idxtype *adjwgt, idxtype *nparts, float *ubvec, 
                   idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_EdgeND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options,
                   idxtype *perm, idxtype *iperm);

void __cdecl METIS_NodeND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options,
                   idxtype *perm, idxtype *iperm);


void __cdecl METIS_NodeWND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *numflag,
                   idxtype *options, idxtype *perm, idxtype *iperm);

void __cdecl METIS_PartGraphKway2(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part);

void __cdecl METIS_WPartGraphKway2(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_NodeNDP(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, idxtype npes, idxtype *options, 
                   idxtype *perm, idxtype *iperm, idxtype *sizes);

void __cdecl METIS_NodeComputeSeparator(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *options, idxtype *sepsize, idxtype *part);

void __cdecl METIS_EdgeComputeSeparator(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *options, idxtype *sepsize, idxtype *part);

void __cdecl METIS_PartGraphRecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part);

void __cdecl METIS_WPartGraphRecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, 
                   idxtype *options, idxtype *edgecut, idxtype *part);

void __cdecl METIS_PartFillGraph(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, 
                   idxtype *edgecut, idxtype *part);


#ifdef __cplusplus
}
#endif




/*------------------------------------------------------------------------
* Constant definitions 
*-------------------------------------------------------------------------*/

/* Matching Schemes */
#define MTYPE_RM		1
#define MTYPE_HEM		2
#define MTYPE_SHEM		3
#define MTYPE_SHEMKWAY		4
#define MTYPE_SHEBM_ONENORM	5
#define MTYPE_SHEBM_INFNORM	6
#define MTYPE_SBHEM_ONENORM	7
#define MTYPE_SBHEM_INFNORM	8

/* Initial partitioning schemes for PMETIS and ONMETIS */
#define ITYPE_GGPKL		1
#define ITYPE_GGPKLNODE		2
#define ITYPE_RANDOM		2

/* Refinement schemes for PMETIS */
#define RTYPE_FM		1

/* Initial partitioning schemes for KMETIS */
#define ITYPE_PMETIS		1

/* Refinement schemes for KMETIS */
#define RTYPE_KWAYRANDOM	1
#define RTYPE_KWAYGREEDY	2
#define RTYPE_KWAYRANDOM_MCONN	3

/* Refinement schemes for ONMETIS */
#define RTYPE_SEP2SIDED		1
#define RTYPE_SEP1SIDED		2

/* Initial Partitioning Schemes for McKMETIS */
#define ITYPE_McPMETIS		1   	/* Simple McPMETIS */
#define ITYPE_McHPMETIS		2	/* horizontally relaxed McPMETIS */


/* Debug Levels */
#define DBG_TIME	1		/* Perform timing analysis */
#define DBG_OUTPUT	2
#define DBG_COARSEN   	4		/* Show the coarsening progress */
#define DBG_REFINE	8		/* Show info on communication during folding */
#define DBG_IPART	16		/* Show info on initial partition */
#define DBG_MOVEINFO	32		/* Show info on communication during folding */
#define DBG_KWAYPINFO	64		/* Show info on communication during folding */
#define DBG_SEPINFO	128		/* Show info on communication during folding */


/* Metis's version number */
#define METIS_VER_MAJOR         5
#define METIS_VER_MINOR         0
#define METIS_VER_SUBMINOR      0



#endif  /* METIS_H */
