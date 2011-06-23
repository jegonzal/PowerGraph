/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef METIS_HPP
#define METIS_HPP

namespace metis {
extern "C" {

void  METIS_EstimateMemory(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *numflag, 
                   int64_t *optype, int64_t *nbytes);

void  METIS_PartGraphKway(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt, 
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, int64_t *options, 
                   int64_t *edgecut, int64_t *part); 

void  METIS_WPartGraphKway(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt, 
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, float *tpwgts, 
                   int64_t *options, int64_t *edgecut, int64_t *part); 

void  METIS_PartGraphVKway(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *vsize, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, int64_t *options, 
                   int64_t *volume, int64_t *part);

void  METIS_WPartGraphVKway(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *vsize, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, float *tpwgts, 
                   int64_t *options, int64_t *volume, int64_t *part);

int64_t   METIS_MeshToDualCount(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *elms, int64_t *etype,
                   int64_t *numflag);

void  METIS_MeshToDual(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *elms, int64_t *etype, 
                   int64_t *numflag,  int64_t *dxadj, int64_t *dadjncy);

int64_t   METIS_MixedMeshToDualCount(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t * elms, 
                   int64_t *etype, int64_t *numflag, int64_t *conmat, int64_t custom);

void  METIS_MixedMeshToDual(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *elms, 
                   int64_t *etype, int64_t *numflag,int64_t *dxadj, int64_t *dadjncy,int64_t *conmat,
                   int64_t custom);

void  METIS_MeshToNodal(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype, int64_t *numflag, 
                   int64_t *dxadj, int64_t *dadjncy);

void  METIS_MixedMeshToNodal(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype,
                   int64_t *numflag, int64_t *dxadj, int64_t *dadjncy);

void  METIS_PartMeshNodal(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype, int64_t *numflag,
                   int64_t *nparts, int64_t *edgecut, int64_t *epart, int64_t *npart);

void  METIS_PartMixedMeshNodal(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype, int64_t *numflag,
                   int64_t *nparts, int64_t *edgecut, int64_t *epart, int64_t *npart);

void  METIS_PartMeshDual(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype, int64_t *numflag, 
                   int64_t *nparts, int64_t *edgecut, int64_t *epart, int64_t *npart, int64_t wgtflag, 
                   int64_t * vwgt);

void  METIS_PartMixedMeshDual(int64_t *ne, int64_t *nn, int64_t *elmnts, int64_t *etype, int64_t *numflag,  
                   int64_t *nparts, int64_t *edgecut, int64_t *epart, int64_t *npart, int64_t *conmat, 
                   int64_t custom, int64_t wgtflag, int64_t *vwgt);

void  METIS_mCPartGraphKway(int64_t *nvtxs, int64_t *ncon, int64_t *xadj, int64_t *adjncy, 
                   int64_t *vwgt, int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, 
                   float *rubvec, int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_mCPartGraphRecursive(int64_t *nvtxs, int64_t *ncon, int64_t *xadj, int64_t *adjncy,
                   int64_t *vwgt, int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts,
                   int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_mCHPartGraphRecursive(int64_t *nvtxs, int64_t *ncon, int64_t *xadj, int64_t *adjncy,
                   int64_t *vwgt, int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts,
                   float *ubvec, int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_mCPartGraphRecursiveInternal(int64_t *nvtxs, int64_t *ncon, int64_t *xadj, 
                   int64_t *adjncy, float *nvwgt, int64_t *adjwgt, int64_t *nparts, int64_t *options, 
                   int64_t *edgecut, int64_t *part);

void  METIS_mCHPartGraphRecursiveInternal(int64_t *nvtxs, int64_t *ncon, int64_t *xadj, 
                   int64_t *adjncy, float *nvwgt, int64_t *adjwgt, int64_t *nparts, float *ubvec, 
                   int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_EdgeND(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *numflag, int64_t *options,
                   int64_t *perm, int64_t *iperm);

void  METIS_NodeND(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *numflag, int64_t *options,
                   int64_t *perm, int64_t *iperm);


void  METIS_NodeWND(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt, int64_t *numflag,
                   int64_t *options, int64_t *perm, int64_t *iperm);

void  METIS_PartGraphKway2(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, int64_t *options, 
                   int64_t *edgecut, int64_t *part);

void  METIS_WPartGraphKway2(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, float *tpwgts, 
                   int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_NodeNDP(int64_t nvtxs, int64_t *xadj, int64_t *adjncy, int64_t npes, int64_t *options, 
                   int64_t *perm, int64_t *iperm, int64_t *sizes);

void  METIS_NodeComputeSeparator(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *options, int64_t *sepsize, int64_t *part);

void  METIS_EdgeComputeSeparator(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *options, int64_t *sepsize, int64_t *part);

void  METIS_PartGraphRecursive(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, int64_t *options, 
                   int64_t *edgecut, int64_t *part);

void  METIS_WPartGraphRecursive(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, float *tpwgts, 
                   int64_t *options, int64_t *edgecut, int64_t *part);

void  METIS_PartFillGraph(int64_t *nvtxs, int64_t *xadj, int64_t *adjncy, int64_t *vwgt,
                   int64_t *adjwgt, int64_t *wgtflag, int64_t *numflag, int64_t *nparts, int64_t *options, 
                   int64_t *edgecut, int64_t *part);
}

typedef int64_t idxtype;


}

#endif

