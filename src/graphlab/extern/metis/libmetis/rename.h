/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains header files
 *
 * Started 10/2/97
 * George
 *
 * $Id: rename.h,v 1.1 2002/08/10 04:34:09 karypis Exp $
 *
 */

/*#define _RENAME_H_*/

#ifndef _RENAME_H_
#define _RENAME_H_

/* balance.c */
#define Balance2Way			libmetis__Balance2Way
#define Bnd2WayBalance			libmetis__Bnd2WayBalance
#define General2WayBalance		libmetis__General2WayBalance


/* bucketsort.c */
#define BucketSortKeysInc		libmetis__BucketSortKeysInc


/* ccgraph.c */
#define CreateCoarseGraph		libmetis__CreateCoarseGraph
#define CreateCoarseGraphNoMask		libmetis__CreateCoarseGraphNoMask
#define CreateCoarseGraph_NVW 		libmetis__CreateCoarseGraph_NVW
#define SetUpCoarseGraph		libmetis__SetUpCoarseGraph
#define ReAdjustMemory			libmetis__ReAdjustMemory


/* coarsen.c */
#define Coarsen2Way			libmetis__Coarsen2Way


/* compress.c */
#define CompressGraph			libmetis__CompressGraph
#define PruneGraph			libmetis__PruneGraph


/* debug.c */
#define ComputeCut			libmetis__ComputeCut
#define CheckBnd			libmetis__CheckBnd
#define CheckBnd2			libmetis__CheckBnd2
#define CheckNodeBnd			libmetis__CheckNodeBnd
#define CheckRInfo			libmetis__CheckRInfo
#define CheckNodePartitionParams	libmetis__CheckNodePartitionParams
#define IsSeparable			libmetis__IsSeparable


/* estmem.c */
#define EstimateCFraction		libmetis__EstimateCFraction
#define ComputeCoarseGraphSize		libmetis__ComputeCoarseGraphSize


/* fm.c */
#define FM_2WayEdgeRefine		libmetis__FM_2WayEdgeRefine


/* fortran.c */
#define Change2CNumbering		libmetis__Change2CNumbering
#define Change2FNumbering		libmetis__Change2FNumbering
#define Change2FNumbering2		libmetis__Change2FNumbering2
#define Change2FNumberingOrder		libmetis__Change2FNumberingOrder
#define ChangeMesh2CNumbering		libmetis__ChangeMesh2CNumbering
#define ChangeMesh2FNumbering		libmetis__ChangeMesh2FNumbering
#define ChangeMesh2FNumbering2		libmetis__ChangeMesh2FNumbering2


/* graph.c */
#define SetUpGraph			libmetis__SetUpGraph
#define SetUpGraph2			libmetis__SetUpGraph2
#define VolSetUpGraph			libmetis__VolSetUpGraph
#define RandomizeGraph			libmetis__RandomizeGraph
#define IsConnectedSubdomain		libmetis__IsConnectedSubdomain
#define IsConnected			libmetis__IsConnected
#define IsConnected2			libmetis__IsConnected2
#define FindComponents			libmetis__FindComponents


/* initpart.c */
#define Init2WayPartition		libmetis__Init2WayPartition
#define InitSeparator			libmetis__InitSeparator
#define GrowBisection			libmetis__GrowBisection
#define GrowBisectionNode		libmetis__GrowBisectionNode
#define RandomBisection			libmetis__RandomBisection


/* kmetis.c */
#define MlevelKWayPartitioning		libmetis__MlevelKWayPartitioning


/* kvmetis.c */
#define MlevelVolKWayPartitioning	libmetis__MlevelVolKWayPartitioning


/* kwayfm.c */
#define Random_KWayEdgeRefine		libmetis__Random_KWayEdgeRefine
#define Greedy_KWayEdgeRefine		libmetis__Greedy_KWayEdgeRefine
#define Greedy_KWayEdgeBalance		libmetis__Greedy_KWayEdgeBalance


/* kwayrefine.c */
#define RefineKWay			libmetis__RefineKWay
#define AllocateKWayPartitionMemory	libmetis__AllocateKWayPartitionMemory
#define ComputeKWayPartitionParams	libmetis__ComputeKWayPartitionParams
#define ProjectKWayPartition		libmetis__ProjectKWayPartition
#define IsBalanced			libmetis__IsBalanced
#define ComputeKWayBoundary		libmetis__ComputeKWayBoundary
#define ComputeKWayBalanceBoundary	libmetis__ComputeKWayBalanceBoundary


/* kwayvolfm.c */
#define Random_KWayVolRefine		libmetis__Random_KWayVolRefine
#define Random_KWayVolRefineMConn	libmetis__Random_KWayVolRefineMConn
#define Greedy_KWayVolBalance		libmetis__Greedy_KWayVolBalance
#define Greedy_KWayVolBalanceMConn	libmetis__Greedy_KWayVolBalanceMConn
#define KWayVolUpdate			libmetis__KWayVolUpdate
#define ComputeKWayVolume		libmetis__ComputeKWayVolume
#define ComputeVolume			libmetis__ComputeVolume
#define CheckVolKWayPartitionParams	libmetis__CheckVolKWayPartitionParams
#define ComputeVolSubDomainGraph	libmetis__ComputeVolSubDomainGraph
#define EliminateVolSubDomainEdges	libmetis__EliminateVolSubDomainEdges


/* kwayvolrefine.c */
#define RefineVolKWay			libmetis__RefineVolKWay
#define AllocateVolKWayPartitionMemory	libmetis__AllocateVolKWayPartitionMemory
#define ComputeVolKWayPartitionParams	libmetis__ComputeVolKWayPartitionParams
#define ComputeKWayVolGains		libmetis__ComputeKWayVolGains
#define ProjectVolKWayPartition		libmetis__ProjectVolKWayPartition
#define ComputeVolKWayBoundary		libmetis__ComputeVolKWayBoundary
#define ComputeVolKWayBalanceBoundary	libmetis__ComputeVolKWayBalanceBoundary


/* match.c */
#define Match_RM			libmetis__Match_RM
#define Match_RM_NVW			libmetis__Match_RM_NVW
#define Match_HEM			libmetis__Match_HEM
#define Match_SHEM			libmetis__Match_SHEM


/* mbalance.c */
#define MocBalance2Way			libmetis__MocBalance2Way
#define MocGeneral2WayBalance		libmetis__MocGeneral2WayBalance


/* mbalance2.c */
#define MocBalance2Way2			libmetis__MocBalance2Way2
#define MocGeneral2WayBalance2		libmetis__MocGeneral2WayBalance2
#define SelectQueue3			libmetis__SelectQueue3


/* mcoarsen.c */
#define MCCoarsen2Way			libmetis__MCCoarsen2Way


/* memory.c */
#define AllocateWorkSpace		libmetis__AllocateWorkSpace
#define FreeWorkSpace			libmetis__FreeWorkSpace
#define WspaceAvail			libmetis__WspaceAvail
#define idxwspacemalloc			libmetis__idxwspacemalloc
#define idxwspacefree			libmetis__idxwspacefree
#define fwspacemalloc			libmetis__fwspacemalloc
#define CreateGraph			libmetis__CreateGraph
#define InitGraph			libmetis__InitGraph
#define FreeRData			libmetis__FreeRData
#define FreeGraph			libmetis__FreeGraph


/* mesh.c */
#define TRIDUALMETIS			libmetis__TRIDUALMETIS
#define TETDUALMETIS			libmetis__TETDUALMETIS
#define HEXDUALMETIS			libmetis__HEXDUALMETIS
#define TRINODALMETIS			libmetis__TRINODALMETIS
#define TETNODALMETIS			libmetis__TETNODALMETIS
#define HEXNODALMETIS			libmetis__HEXNODALMETIS


/* mfm.c */
#define MocFM_2WayEdgeRefine		libmetis__MocFM_2WayEdgeRefine
#define SelectQueue			libmetis__SelectQueue
#define BetterBalance			libmetis__BetterBalance
#define Compute2WayHLoadImbalance	libmetis__Compute2WayHLoadImbalance
#define Compute2WayHLoadImbalanceVec	libmetis__Compute2WayHLoadImbalanceVec


/* mfm2.c */
#define MocFM_2WayEdgeRefine2		libmetis__MocFM_2WayEdgeRefine2
#define SelectQueue2			libmetis__SelectQueue2
#define IsBetter2wayBalance		libmetis__IsBetter2wayBalance


/* mincover.c */
#define MinCover			libmetis__MinCover
#define MinCover_Augment		libmetis__MinCover_Augment
#define MinCover_Decompose		libmetis__MinCover_Decompose
#define MinCover_ColDFS			libmetis__MinCover_ColDFS
#define MinCover_RowDFS			libmetis__MinCover_RowDFS


/* minitpart.c */
#define MocInit2WayPartition		libmetis__MocInit2WayPartition
#define MocGrowBisection		libmetis__MocGrowBisection
#define MocRandomBisection		libmetis__MocRandomBisection
#define MocInit2WayBalance		libmetis__MocInit2WayBalance
#define SelectQueueoneWay		libmetis__SelectQueueoneWay


/* minitpart2.c */
#define MocInit2WayPartition2		libmetis__MocInit2WayPartition2
#define MocGrowBisection2		libmetis__MocGrowBisection2
#define MocGrowBisectionNew2		libmetis__MocGrowBisectionNew2
#define MocInit2WayBalance2		libmetis__MocInit2WayBalance2
#define SelectQueueOneWay2		libmetis__SelectQueueOneWay2


/* mkmetis.c */
#define MCMlevelKWayPartitioning	libmetis__MCMlevelKWayPartitioning


/* mkwayfmh.c */
#define MCRandom_KWayEdgeRefineHorizontal	libmetis__MCRandom_KWayEdgeRefineHorizontal
#define MCGreedy_KWayEdgeBalanceHorizontal	libmetis__MCGreedy_KWayEdgeBalanceHorizontal
#define AreAllHVwgtsBelow			libmetis__AreAllHVwgtsBelow
#define AreAllHVwgtsAbove			libmetis__AreAllHVwgtsAbove
#define ComputeHKWayLoadImbalance		libmetis__ComputeHKWayLoadImbalance
#define MocIsHBalanced				libmetis__MocIsHBalanced
#define IsHBalanceBetterFT			libmetis__IsHBalanceBetterFT
#define IsHBalanceBetterTT			libmetis__IsHBalanceBetterTT


/* mkwayrefine.c */
#define MocRefineKWayHorizontal		libmetis__MocRefineKWayHorizontal
#define MocAllocateKWayPartitionMemory	libmetis__MocAllocateKWayPartitionMemory
#define MocComputeKWayPartitionParams	libmetis__MocComputeKWayPartitionParams
#define MocProjectKWayPartition		libmetis__MocProjectKWayPartition
#define MocComputeKWayBalanceBoundary	libmetis__MocComputeKWayBalanceBoundary


/* mmatch.c */
#define MCMatch_RM			libmetis__MCMatch_RM
#define MCMatch_HEM			libmetis__MCMatch_HEM
#define MCMatch_SHEM			libmetis__MCMatch_SHEM
#define MCMatch_SHEBM			libmetis__MCMatch_SHEBM
#define MCMatch_SBHEM			libmetis__MCMatch_SBHEM
#define BetterVBalance			libmetis__BetterVBalance
#define AreAllVwgtsBelowFast		libmetis__AreAllVwgtsBelowFast


/* mmd.c */
#define genmmd				libmetis__genmmd
#define mmdelm				libmetis__mmdelm
#define mmdint				libmetis__mmdint
#define mmdnum				libmetis__mmdnum
#define mmdupd				libmetis__mmdupd


/* mpmetis.c */
#define MCMlevelRecursiveBisection	libmetis__MCMlevelRecursiveBisection
#define MCHMlevelRecursiveBisection	libmetis__MCHMlevelRecursiveBisection
#define MCMlevelEdgeBisection		libmetis__MCMlevelEdgeBisection
#define MCHMlevelEdgeBisection		libmetis__MCHMlevelEdgeBisection


/* mrefine.c */
#define MocRefine2Way			libmetis__MocRefine2Way
#define MocAllocate2WayPartitionMemory	libmetis__MocAllocate2WayPartitionMemory
#define MocCompute2WayPartitionParams	libmetis__MocCompute2WayPartitionParams
#define MocProject2WayPartition		libmetis__MocProject2WayPartition


/* mrefine2.c */
#define MocRefine2Way2			libmetis__MocRefine2Way2


/* mutil.c */
#define AreAllVwgtsBelow		libmetis__AreAllVwgtsBelow
#define AreAnyVwgtsBelow		libmetis__AreAnyVwgtsBelow
#define AreAllVwgtsAbove		libmetis__AreAllVwgtsAbove
#define ComputeLoadImbalance		libmetis__ComputeLoadImbalance
#define AreAllBelow			libmetis__AreAllBelow


/* myqsort.c */
#define iidxsort			libmetis__iidxsort
#define ikeysort			libmetis__ikeysort
#define ikeyvalsort			libmetis__ikeyvalsort
#define idkeysort			libmetis__idkeysort


/* ometis.c */
#define MlevelNestedDissection		libmetis__MlevelNestedDissection
#define MlevelNestedDissectionCC	libmetis__MlevelNestedDissectionCC
#define MlevelNodeBisectionMultiple	libmetis__MlevelNodeBisectionMultiple
#define MlevelNodeBisection		libmetis__MlevelNodeBisection
#define SplitGraphOrder			libmetis__SplitGraphOrder
#define MMDOrder			libmetis__MMDOrder
#define SplitGraphOrderCC		libmetis__SplitGraphOrderCC


/* parmetis.c */
#define MlevelNestedDissectionP		libmetis__MlevelNestedDissectionP


/* pmetis.c */
#define MlevelRecursiveBisection	libmetis__MlevelRecursiveBisection
#define MlevelEdgeBisection		libmetis__MlevelEdgeBisection
#define SplitGraphPart			libmetis__SplitGraphPart
#define SetUpSplitGraph			libmetis__SetUpSplitGraph


/* pqueue.c */
#define PQueueInit			libmetis__PQueueInit
#define PQueueReset			libmetis__PQueueReset
#define PQueueFree			libmetis__PQueueFree
#define PQueueInsert			libmetis__PQueueInsert
#define PQueueDelete			libmetis__PQueueDelete
#define PQueueUpdate			libmetis__PQueueUpdate
#define PQueueUpdateUp			libmetis__PQueueUpdateUp
#define PQueueGetMax			libmetis__PQueueGetMax
#define PQueueSeeMax			libmetis__PQueueSeeMax
#define CheckHeap			libmetis__CheckHeap


/* refine.c */
#define Refine2Way			libmetis__Refine2Way
#define Allocate2WayPartitionMemory	libmetis__Allocate2WayPartitionMemory
#define Compute2WayPartitionParams	libmetis__Compute2WayPartitionParams
#define Project2WayPartition		libmetis__Project2WayPartition


/* separator.c */
#define ConstructSeparator		libmetis__ConstructSeparator
#define ConstructMinCoverSeparator0	libmetis__ConstructMinCoverSeparator0
#define ConstructMinCoverSeparator	libmetis__ConstructMinCoverSeparator


/* sfm.c */
#define FM_2WayNodeRefine		libmetis__FM_2WayNodeRefine
#define FM_2WayNodeRefineEqWgt		libmetis__FM_2WayNodeRefineEqWgt
#define FM_2WayNodeRefine_OneSided	libmetis__FM_2WayNodeRefine_OneSided
#define FM_2WayNodeBalance		libmetis__FM_2WayNodeBalance
#define ComputeMaxNodeGain		libmetis__ComputeMaxNodeGain


/* srefine.c */
#define Refine2WayNode			libmetis__Refine2WayNode
#define Allocate2WayNodePartitionMemory	libmetis__Allocate2WayNodePartitionMemory
#define Compute2WayNodePartitionParams	libmetis__Compute2WayNodePartitionParams
#define Project2WayNodePartition	libmetis__Project2WayNodePartition


/* stat.c */
#define ComputePartitionInfo		libmetis__ComputePartitionInfo
#define ComputePartitionBalance		libmetis__ComputePartitionBalance
#define ComputeElementBalance		libmetis__ComputeElementBalance


/* subdomains.c */
#define Random_KWayEdgeRefineMConn	libmetis__Random_KWayEdgeRefineMConn
#define Greedy_KWayEdgeBalanceMConn	libmetis__Greedy_KWayEdgeBalanceMConn
#define PrintSubDomainGraph		libmetis__PrintSubDomainGraph
#define ComputeSubDomainGraph		libmetis__ComputeSubDomainGraph
#define EliminateSubDomainEdges		libmetis__EliminateSubDomainEdges
#define MoveGroupMConn			libmetis__MoveGroupMConn
#define EliminateComponents		libmetis__EliminateComponents
#define MoveGroup			libmetis__MoveGroup


/* timing.c */
#define InitTimers			libmetis__InitTimers
#define PrintTimers			libmetis__PrintTimers


/* util.c */
#define idxmalloc			libmetis__idxmalloc
#define idxsmalloc			libmetis__idxsmalloc
#define idxset				libmetis__idxset
#define idxargmax		        libmetis__idxargmax
#define idxargmin			libmetis__idxargmin
#define idxsum				libmetis__idxsum
#define idxaxpy				libmetis__idxaxpy
#define idxargmax_strd			libmetis__idxargmax_strd
#define famax2				libmetis__famax2
#define RandomPermute			libmetis__RandomPermute
#define InitRandom			libmetis__InitRandom


#endif


