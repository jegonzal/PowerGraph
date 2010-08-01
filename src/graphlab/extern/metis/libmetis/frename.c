/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * frename.c
 * 
 * This file contains some renaming routines to deal with different Fortran compilers
 *
 * Started 9/15/97
 * George
 *
 */

#include <metislib.h>


void METIS_PARTGRAPHRECURSIVE(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{ 
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part); 
}
void metis_partgraphrecursive_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_WPARTGRAPHRECURSIVE(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_PARTGRAPHKWAY(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}



void METIS_WPARTGRAPHKWAY(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_EDGEND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEWND(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *numflag, idxtype *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}


#ifdef XXXX
void METIS_PARTMESHNODAL(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal_(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal__(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}


void METIS_PARTMESHDUAL(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart, 0, NULL);
}
void metis_partmeshdual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart, 0, NULL);
}
void metis_partmeshdual_(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart, 0, NULL);
}
void metis_partmeshdual__(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *nparts, idxtype *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart, 0, NULL);
}


void METIS_MESHTONODAL(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal_(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal__(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}


void METIS_MESHTODUAL(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual_(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual__(idxtype *ne, idxtype *nn, idxtype *elmnts, idxtype *etype, idxtype *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
#endif

void METIS_ESTIMATEMEMORY(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *optype, idxtype *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *optype, idxtype *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *optype, idxtype *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *numflag, idxtype *optype, idxtype *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}



void METIS_MCPARTGRAPHRECURSIVE(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive_(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive__(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_MCPARTGRAPHKWAY(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway_(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway__(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}


void METIS_PARTGRAPHVKWAY(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkaway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkaway_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkaway__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}

void METIS_WPARTGRAPHVKWAY(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkaway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkaway_(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkaway__(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}



