/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.27 2003/05/03 16:10:48 karypis Exp $
 *
 */

#ifndef _LIBMETIS_PROTO_H_
#define _LIBMETIS_PROTO_H_

/* balance.c */
void Balance2Way(CtrlType *, GraphType *, idxtype *, float);
void Bnd2WayBalance(CtrlType *, GraphType *, idxtype *);
void General2WayBalance(CtrlType *, GraphType *, idxtype *);

/* bucketsort.c */
void BucketSortKeysInc(idxtype, idxtype, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void CreateCoarseGraph(CtrlType *, GraphType *, idxtype, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, idxtype, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, idxtype, idxtype *, idxtype *);
GraphType *SetUpCoarseGraph(GraphType *, idxtype, idxtype);
void ReAdjustMemory(GraphType *, GraphType *, idxtype);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* compress.c */
void CompressGraph(CtrlType *, GraphType *, idxtype, idxtype *, idxtype *, idxtype *, idxtype *);
void PruneGraph(CtrlType *, GraphType *, idxtype, idxtype *, idxtype *, idxtype *, float);

/* debug.c */
idxtype ComputeCut(GraphType *, idxtype *);
idxtype ComputeMaxCut(GraphType *graph, idxtype nparts, idxtype *where);
idxtype CheckBnd(GraphType *);
idxtype CheckBnd2(GraphType *);
idxtype CheckNodeBnd(GraphType *, idxtype);
idxtype CheckRInfo(RInfoType *);
idxtype CheckNodePartitionParams(GraphType *);
idxtype IsSeparable(GraphType *);

/* estmem.c */
void EstimateCFraction(idxtype, idxtype *, idxtype *, float *, float *);
idxtype ComputeCoarseGraphSize(idxtype, idxtype *, idxtype *, idxtype, idxtype *, idxtype *, idxtype *);

/* fm.c */
void FM_2WayEdgeRefine(CtrlType *, GraphType *, idxtype *, idxtype);

/* fortran.c */
void Change2CNumbering(idxtype, idxtype *, idxtype *);
void Change2FNumbering(idxtype, idxtype *, idxtype *, idxtype *);
void Change2FNumbering2(idxtype, idxtype *, idxtype *);
void Change2FNumberingOrder(idxtype, idxtype *, idxtype *, idxtype *, idxtype *);
void ChangeMesh2CNumbering(idxtype, idxtype *);
void ChangeMesh2FNumbering(idxtype, idxtype *, idxtype, idxtype *, idxtype *);
void ChangeMesh2FNumbering2(idxtype, idxtype *, idxtype, idxtype, idxtype *, idxtype *);
void ChangeMesh2FNumbering3(idxtype, idxtype *);

/* frename.c */
void METIS_PARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void METIS_PARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_partgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_WPARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void metis_wpartgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *); 
void METIS_EDGEND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_edgend__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_NODEND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodend__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_NODEWND(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void metis_nodewnd__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partmeshnodal__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *);
void metis_partmeshdual(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *);
void metis_partmeshdual_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *);
void metis_partmeshdual__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *);
void METIS_MESHTONODAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtonodal__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MESHTODUAL(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_meshtodual__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_estimatememory__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MCPARTGRAPHRECURSIVE(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphrecursive__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_MCPARTGRAPHKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_mcpartgraphkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void METIS_PARTGRAPHVKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void metis_partgraphvkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void METIS_WPARTGRAPHVKWAY(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway_(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);
void metis_wpartgraphvkway__(idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, float *, idxtype *, idxtype *, idxtype *);

/* graph.c */
void SetUpGraph(GraphType *, idxtype, idxtype, idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype);
void SetUpGraph2(GraphType *, idxtype, idxtype, idxtype *, idxtype *, float *, idxtype *);
void VolSetUpGraph(GraphType *, idxtype, idxtype, idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype);
void RandomizeGraph(GraphType *);
idxtype IsConnectedSubdomain(CtrlType *, GraphType *, idxtype, idxtype);
idxtype IsConnected(CtrlType *, GraphType *, idxtype);
idxtype IsConnected2(GraphType *, idxtype);
idxtype FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);

/* initpart.c */
void Init2WayPartition(CtrlType *, GraphType *, idxtype *, float);
void InitSeparator(CtrlType *, GraphType *, float);
void GrowBisection(CtrlType *, GraphType *, idxtype *, float);
void GrowBisectionNode(CtrlType *, GraphType *, float);
void RandomBisection(CtrlType *, GraphType *, idxtype *, float);

/* kmetis.c */
idxtype MlevelKWayPartitioning(CtrlType *, GraphType *, idxtype, idxtype *, float *, float);

/* kvmetis.c */
idxtype MlevelVolKWayPartitioning(CtrlType *, GraphType *, idxtype, idxtype *, float *, float);

/* kwayfm.c */
void Random_KWayEdgeRefine(CtrlType *, GraphType *, idxtype, float *, float, idxtype, idxtype);
void Greedy_KWayEdgeRefine(CtrlType *, GraphType *, idxtype, float *, float, idxtype);
void Greedy_KWayEdgeBalance(CtrlType *, GraphType *, idxtype, float *, float, idxtype);

/* kwayrefine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, idxtype, float *, float);
void RefineKWayRefinement(CtrlType *, GraphType *, GraphType *, idxtype, float *, float);
void AllocateKWayPartitionMemory(CtrlType *, GraphType *, idxtype);
void ComputeKWayPartitionParams(CtrlType *, GraphType *, idxtype);
void ProjectKWayPartition(CtrlType *, GraphType *, idxtype);
idxtype IsBalanced(idxtype *, idxtype, float *, float);
void ComputeKWayBoundary(CtrlType *, GraphType *, idxtype);
void ComputeKWayBalanceBoundary(CtrlType *, GraphType *, idxtype);

/* kwayvolfm.c */
void Random_KWayVolRefine(CtrlType *, GraphType *, idxtype, float *, float, idxtype, idxtype);
void Random_KWayVolRefineMConn(CtrlType *, GraphType *, idxtype, float *, float, idxtype, idxtype);
void Greedy_KWayVolBalance(CtrlType *, GraphType *, idxtype, float *, float, idxtype);
void Greedy_KWayVolBalanceMConn(CtrlType *, GraphType *, idxtype, float *, float, idxtype);
void KWayVolUpdate(CtrlType *, GraphType *, idxtype, idxtype, idxtype, idxtype *, idxtype *, idxtype *);
void ComputeKWayVolume(GraphType *, idxtype, idxtype *, idxtype *, idxtype *);
idxtype ComputeVolume(GraphType *, idxtype *);
void CheckVolKWayPartitionParams(CtrlType *, GraphType *, idxtype);
void ComputeVolSubDomainGraph(GraphType *, idxtype, idxtype *, idxtype *);
void EliminateVolSubDomainEdges(CtrlType *, GraphType *, idxtype, float *);
void EliminateVolComponents(CtrlType *, GraphType *, idxtype, float *, float);

/* kwayvolrefine.c */
void RefineVolKWay(CtrlType *, GraphType *, GraphType *, idxtype, float *, float);
void AllocateVolKWayPartitionMemory(CtrlType *, GraphType *, idxtype);
void ComputeVolKWayPartitionParams(CtrlType *, GraphType *, idxtype);
void ComputeKWayVolGains(CtrlType *, GraphType *, idxtype);
void ProjectVolKWayPartition(CtrlType *, GraphType *, idxtype);
void ComputeVolKWayBoundary(CtrlType *, GraphType *, idxtype);
void ComputeVolKWayBalanceBoundary(CtrlType *, GraphType *, idxtype);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);

/* mbalance.c */
void MocBalance2Way(CtrlType *, GraphType *, float *, float);
void MocGeneral2WayBalance(CtrlType *, GraphType *, float *, float);

/* mbalance2.c */
void MocBalance2Way2(CtrlType *, GraphType *, float *, float *);
void MocGeneral2WayBalance2(CtrlType *, GraphType *, float *, float *);
void SelectQueue3(idxtype, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2], float *);

/* mcoarsen.c */
GraphType *MCCoarsen2Way(CtrlType *, GraphType *);

/* memory.c */
void AllocateWorkSpace(CtrlType *, GraphType *, idxtype);
void FreeWorkSpace(CtrlType *, GraphType *);
idxtype WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, idxtype);
void idxwspacefree(CtrlType *, idxtype);
float *fwspacemalloc(CtrlType *, idxtype);
void fwspacefree(CtrlType *, idxtype);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeRData(GraphType *);
void FreeGraph(GraphType *, int flag);

/* mesh.c */
idxtype GENDUALMETIS_COUNT(idxtype , idxtype , idxtype , idxtype *, idxtype *);
void GENDUALMETIS(idxtype, idxtype, idxtype, idxtype *, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(idxtype, idxtype, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(idxtype, idxtype, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(idxtype, idxtype, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(idxtype, idxtype, idxtype *, idxtype *, idxtype *adjncy);
void LINENODALMETIS(idxtype, idxtype, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */

/* mfm.c */
void MocFM_2WayEdgeRefine(CtrlType *, GraphType *, float *, idxtype);
void SelectQueue(idxtype, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2]);
idxtype BetterBalance(idxtype, float *, float *, float *);
float Compute2WayHLoadImbalance(idxtype, float *, float *);
void Compute2WayHLoadImbalanceVec(idxtype, float *, float *, float *);

/* mfm2.c */
void MocFM_2WayEdgeRefine2(CtrlType *, GraphType *, float *, float *, idxtype);
void SelectQueue2(idxtype, float *, float *, idxtype *, idxtype *, PQueueType [MAXNCON][2], float *);
idxtype IsBetter2wayBalance(idxtype, float *, float *, float *);

/* mincover.o */
void MinCover(idxtype *, idxtype *, idxtype, idxtype, idxtype *, idxtype *);
idxtype MinCover_Augment(idxtype *, idxtype *, idxtype, idxtype *, idxtype *, idxtype *, idxtype);
void MinCover_Decompose(idxtype *, idxtype *, idxtype, idxtype, idxtype *, idxtype *, idxtype *);
void MinCover_ColDFS(idxtype *, idxtype *, idxtype, idxtype *, idxtype *, idxtype);
void MinCover_RowDFS(idxtype *, idxtype *, idxtype, idxtype *, idxtype *, idxtype);

/* minitpart.c */
void MocInit2WayPartition(CtrlType *, GraphType *, float *, float);
void MocGrowBisection(CtrlType *, GraphType *, float *, float);
void MocRandomBisection(CtrlType *, GraphType *, float *, float);
void MocInit2WayBalance(CtrlType *, GraphType *, float *);
idxtype SelectQueueOneWay(idxtype ncon, float *npwgts, float *tpwgts, idxtype from, PQueueType queues[MAXNCON][2]);

/* minitpart2.c */
void MocInit2WayPartition2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisection2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisectionNew2(CtrlType *, GraphType *, float *, float *);
void MocInit2WayBalance2(CtrlType *, GraphType *, float *, float *);
idxtype SelectQueueOneWay2(idxtype ncon, float *pto, PQueueType queues[MAXNCON][2], float *ubvec);



/* mkmetis.c */
idxtype MCMlevelKWayPartitioning(CtrlType *, GraphType *, idxtype, idxtype *, float *);

/* mkwayfmh.c */
void MCRandom_KWayEdgeRefineHorizontal(CtrlType *, GraphType *, idxtype, float *, idxtype);
void MCGreedy_KWayEdgeBalanceHorizontal(CtrlType *, GraphType *, idxtype, float *, idxtype);
idxtype AreAllHVwgtsBelow(idxtype, float, float *, float, float *, float *);
idxtype AreAllHVwgtsAbove(idxtype, float, float *, float, float *, float *);
void ComputeHKWayLoadImbalance(idxtype, idxtype, float *, float *);
idxtype MocIsHBalanced(idxtype, idxtype, float *, float *);
idxtype IsHBalanceBetterFT(idxtype, idxtype, float *, float *, float *, float *);
idxtype IsHBalanceBetterTT(idxtype, idxtype, float *, float *, float *, float *);

/* mkwayrefine.c */
void MocRefineKWayHorizontal(CtrlType *, GraphType *, GraphType *, idxtype, float *);
void MocAllocateKWayPartitionMemory(CtrlType *, GraphType *, idxtype);
void MocComputeKWayPartitionParams(CtrlType *, GraphType *, idxtype);
void MocProjectKWayPartition(CtrlType *, GraphType *, idxtype);
void MocComputeKWayBalanceBoundary(CtrlType *, GraphType *, idxtype);

/* mmatch.c */
void MCMatch_RM(CtrlType *, GraphType *);
void MCMatch_HEM(CtrlType *, GraphType *);
void MCMatch_SHEM(CtrlType *, GraphType *);
void MCMatch_SHEBM(CtrlType *, GraphType *, idxtype);
void MCMatch_SBHEM(CtrlType *, GraphType *, idxtype);
float BetterVBalance(idxtype, idxtype, float *, float *, float *);
idxtype AreAllVwgtsBelowFast(idxtype, float *, float *, float);

/* mmd.c */
void genmmd(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype , idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *);
void mmdelm(idxtype, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype);
idxtype mmdint(idxtype, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(idxtype, idxtype *, idxtype *, idxtype *);
void mmdupd(idxtype, idxtype, idxtype *, idxtype *, idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype *tag);

/* mpmetis.c */
idxtype MCMlevelRecursiveBisection(CtrlType *, GraphType *, idxtype, idxtype *, float, idxtype);
idxtype MCHMlevelRecursiveBisection(CtrlType *, GraphType *, idxtype, idxtype *, float *, idxtype);
void MCMlevelEdgeBisection(CtrlType *, GraphType *, float *, float);
void MCHMlevelEdgeBisection(CtrlType *, GraphType *, float *, float *);

/* mrefine.c */
void MocRefine2Way(CtrlType *, GraphType *, GraphType *, float *, float);
void MocAllocate2WayPartitionMemory(CtrlType *, GraphType *);
void MocCompute2WayPartitionParams(CtrlType *, GraphType *);
void MocProject2WayPartition(CtrlType *, GraphType *);

/* mrefine2.c */
void MocRefine2Way2(CtrlType *, GraphType *, GraphType *, float *, float *);

/* mutil.c */
idxtype AreAllVwgtsBelow(idxtype, float, float *, float, float *, float);
idxtype AreAnyVwgtsBelow(idxtype, float, float *, float, float *, float);
idxtype AreAllVwgtsAbove(idxtype, float, float *, float, float *, float);
float ComputeLoadImbalance(idxtype, idxtype, float *, float *);
idxtype AreAllBelow(idxtype, float *, float *);

/* myqsort.c */
void iidxsort(size_t, idxtype *);
void ikeysort(size_t, KeyValueType *);
void ikeyvalsort(size_t, KeyValueType *);
void idkeysort(size_t, DKeyValueType *);


/* ometis.c */
void MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, idxtype);
void MlevelNestedDissectionCC(CtrlType *, GraphType *, idxtype *, float, idxtype);
void MlevelNodeBisectionMultiple(CtrlType *, GraphType *, idxtype *, float);
void MlevelNodeBisection(CtrlType *, GraphType *, idxtype *, float);
void SplitGraphOrder(CtrlType *, GraphType *, GraphType *, GraphType *);
void MMDOrder(CtrlType *, GraphType *, idxtype *, idxtype);
idxtype SplitGraphOrderCC(CtrlType *, GraphType *, GraphType *, idxtype, idxtype *, idxtype *);

/* parmetis.c */
void MlevelNestedDissectionP(CtrlType *, GraphType *, idxtype *, idxtype, idxtype, idxtype, idxtype *);

/* pmetis.c */
idxtype MlevelRecursiveBisection(CtrlType *, GraphType *, idxtype, idxtype *, float *, float, idxtype);
void MlevelEdgeBisection(CtrlType *, GraphType *, idxtype *, float);
void SplitGraphPart(CtrlType *, GraphType *, GraphType *, GraphType *);
void SetUpSplitGraph(GraphType *, GraphType *, idxtype, idxtype);

/* pqueue.c */
void PQueueInit(CtrlType *ctrl, PQueueType *, idxtype, idxtype);
void PQueueReset(PQueueType *);
void PQueueFree(CtrlType *ctrl, PQueueType *);
idxtype PQueueGetSize(PQueueType *);
idxtype PQueueInsert(PQueueType *, idxtype, idxtype);
idxtype PQueueDelete(PQueueType *, idxtype, idxtype);
idxtype PQueueUpdate(PQueueType *, idxtype, idxtype, idxtype);
void PQueueUpdateUp(PQueueType *, idxtype, idxtype, idxtype);
idxtype PQueueGetMax(PQueueType *);
idxtype PQueueSeeMax(PQueueType *);
idxtype PQueueGetKey(PQueueType *);
idxtype CheckHeap(PQueueType *);

/* refine.c */
void Refine2Way(CtrlType *, GraphType *, GraphType *, idxtype *, float ubfactor);
void Allocate2WayPartitionMemory(CtrlType *, GraphType *);
void Compute2WayPartitionParams(CtrlType *, GraphType *);
void Project2WayPartition(CtrlType *, GraphType *);

/* separator.c */
void ConstructSeparator(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator0(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator(CtrlType *, GraphType *, float);

/* sfm.c */
void FM_2WayNodeRefine(CtrlType *, GraphType *, float, idxtype);
void FM_2WayNodeRefineEqWgt(CtrlType *, GraphType *, idxtype);
void FM_2WayNodeRefine_OneSided(CtrlType *, GraphType *, float, idxtype);
void FM_2WayNodeBalance(CtrlType *, GraphType *, float);
idxtype ComputeMaxNodeGain(idxtype, idxtype *, idxtype *, idxtype *);

/* srefine.c */
void Refine2WayNode(CtrlType *, GraphType *, GraphType *, float);
void Allocate2WayNodePartitionMemory(CtrlType *, GraphType *);
void Compute2WayNodePartitionParams(CtrlType *, GraphType *);
void Project2WayNodePartition(CtrlType *, GraphType *);

/* stat.c */
void ComputePartitionInfo(GraphType *, idxtype, idxtype *);
void ComputePartitionInfoBipartite(GraphType *, idxtype, idxtype *);
void ComputePartitionBalance(GraphType *, idxtype, idxtype *, float *);
float ComputeElementBalance(idxtype, idxtype, idxtype *);

/* subdomains.c */
void Random_KWayEdgeRefineMConn(CtrlType *, GraphType *, idxtype, float *, float, idxtype, idxtype);
void Greedy_KWayEdgeBalanceMConn(CtrlType *, GraphType *, idxtype, float *, float, idxtype);
void PrintSubDomainGraph(GraphType *, idxtype, idxtype *);
void ComputeSubDomainGraph(GraphType *, idxtype, idxtype *, idxtype *);
void EliminateSubDomainEdges(CtrlType *, GraphType *, idxtype, float *);
void MoveGroupMConn(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype, idxtype, idxtype, idxtype *);
void EliminateComponents(CtrlType *, GraphType *, idxtype, float *, float);
void MoveGroup(CtrlType *, GraphType *, idxtype, idxtype, idxtype, idxtype *, idxtype *);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);

/* util.c */
GK_XMALLOC_PROTO(idxmalloc, idxtype)
GK_XREALLOC_PROTO(idxrealloc, idxtype)
GK_XSMALLOC_PROTO(idxsmalloc, idxtype)
GK_SET_PROTO(idxset, idxtype)
GK_ARGMAX_PROTO(idxargmax, idxtype)
GK_ARGMIN_PROTO(idxargmin, idxtype)
GK_SUM_PROTO(idxsum, idxtype, idxtype)
GK_AXPY_PROTO(idxaxpy, idxtype)

idxtype idxargmax_strd(size_t, idxtype *, idxtype);
idxtype famax2(size_t, float *);
void RandomPermute(size_t, idxtype *, idxtype);
void InitRandom(idxtype);
idxtype strtoidx(const char *nptr, char **endptr, int base);












/***************************************************************
* Programs Directory
****************************************************************/

/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
idxtype smbfct(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);

/* kfmetis.c */
void BalanceFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part);
GraphType *ExtractPartitionGraph(GraphType *graph, idxtype *part, idxtype pid, idxtype *vmap, idxtype *vimap);
void ComputePartitionFillIn(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, idxtype *spart, idxtype pid, idxtype *r_fill, idxtype *r_subfill);
void RefineTopLevelSeparators(CtrlType *ctrl, GraphType *graph, idxtype nparts,  idxtype *part, idxtype *spart, idxtype *, idxtype *, idxtype *, idxtype *);


/* rkmetis.c */
void METIS_RefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, idxtype *options, idxtype *edgecut, idxtype *part);
void METIS_WRefineGraphKway(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag, idxtype *nparts, float *tpwgts, idxtype *options, idxtype *edgecut, idxtype *part);
idxtype MlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part, float *tpwgts, float ubfactor);


/* cmetis.c */
void *METIS_PartGraphForContact(idxtype *nvtxs, idxtype *xadj, idxtype *adjncy,
                double *xyzcoords, idxtype *sflag, idxtype *numflag, idxtype *nparts,
                idxtype *options, idxtype *edgecut, idxtype *part);
void METIS_UpdateContactInfo(void *raw_cinfo, idxtype *nvtxs, double *xyzcoords, idxtype *sflag);
void *METIS_SetupContact0(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part);
void *METIS_SetupContact(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part);
void METIS_FindContacts(void *raw_cinfo, idxtype *nboxes, double *boxcoords, idxtype *nparts, 
               idxtype **r_cntptr, idxtype **r_cntind);
void METIS_FreeContactInfo(void *raw_cinfo);
GraphType *CreatePartitionGraphForContact(idxtype nvtxs, idxtype *xadj, idxtype *adjncy, 
                idxtype *vwgt, idxtype *adjwgt, idxtype cnvtxs, idxtype *part);
idxtype InduceDecissionTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype *sflag, idxtype nparts, 
          idxtype *part, idxtype maxnvtxs, idxtype minnvtxs, float minfrac, idxtype *r_nnodes, idxtype *r_nlnodes,
          DTreeNodeType *dtree, idxtype *dtpart, idxtype *dtipart, idxtype *r_nclean,
          idxtype *r_naclean, idxtype *r_ndirty, idxtype *r_maxdepth, idxtype *marker);
idxtype FindBoxContacts(ContactInfoType *cinfo, double *coords, idxtype *stack, idxtype *cntind, idxtype *marker);
void BuildDTLeafContents(ContactInfoType *cinfo, idxtype *sflag);
void  CheckDTree(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo);
void  CheckDTreeSurface(idxtype nvtxs, double *xyzcoords, idxtype *part, ContactInfoType *cinfo, idxtype *sflag);
void *METIS_PartSurfForContactRCB(idxtype *nvtxs, double *xyzcoords, idxtype *sflag,
                idxtype *nparts, idxtype *part, idxtype *bestdims);
idxtype InduceRCBTree(idxtype nvtxs, DKeyValueType **xyzcand, idxtype firstPID, idxtype nparts,
          idxtype *r_nnodes, idxtype *r_nlnodes, DTreeNodeType *dtree, idxtype *leafpart,
          idxtype *part, idxtype *marker, idxtype *oldBestDims);



/* mrkmetis.c */
void METIS_mCRefineGraphKway(idxtype *nvtxs, idxtype *ncon, idxtype *xadj, idxtype *adjncy,
                          idxtype *vwgt, idxtype *adjwgt, idxtype *wgtflag, idxtype *numflag,
                          idxtype *nparts, float *rubvec, idxtype *options, idxtype *edgecut,
                          idxtype *part);
idxtype MCMlevelKWayRefinement(CtrlType *ctrl, GraphType *graph, idxtype nparts, idxtype *part,
      float *rubvec);


/* streamio.c */
int mprintf(char *format,...);
int msprintf(char *str, char *format,...);
int mfprintf(FILE *stream, char *format,...);
int mscanf(char *format,...);
int msscanf(char *str, char *format,...);
int mfscanf(FILE *stream, char *format,...);



#endif
