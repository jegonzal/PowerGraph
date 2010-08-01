/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h,v 1.11 2003/04/23 02:41:59 karypis Exp $
 */


#define MAXIDX	((idxtype)1<<(8*sizeof(idxtype))-2)

/*************************************************************************
* The following data structure stores idxtype-val - double-key pairs
**************************************************************************/
typedef struct {
  idxtype dim;
  double value;
  idxtype nvtxs;
  idxtype nsvtxs;
  idxtype leafid;
  idxtype partid;

  idxtype left, right;
} DTreeNodeType;



/*************************************************************************
* The following data structure stores idxtype-val - double-key pairs
**************************************************************************/
typedef struct {
  idxtype nvtxs;
  idxtype nnodes, nleafs;
  idxtype *leafptr, *leafind, *leafwgt; /* Information as to which partitions are assigned to the leafs */
  idxtype *part;                        /* The actual partitioning vector */
  idxtype *leafpart;                    /* The leaf-based partitioning vector */
  DTreeNodeType *dtree;                 /* The decission tree itself */
} ContactInfoType;



/*************************************************************************
* The following data structure stores idxtype-key - double-value pairs
**************************************************************************/
typedef struct {
  double key;
  idxtype val;
} DKeyValueType;


/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;


/*************************************************************************
* The following data structure will hold a node of a doubly-linked list.
**************************************************************************/
struct ListNodeType {
  idxtype id;                       	/* The id value of the node */
  struct ListNodeType *prev, *next;     /* It's a doubly-linked list */
};

typedef struct ListNodeType ListNodeType;



/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct PQueueType {
  idxtype type;                     /* The type of the representation used */
  idxtype nnodes;
  idxtype maxnodes;
  idxtype mustfree;

  /* Linear array version of the data structures */
  idxtype pgainspan, ngainspan;     /* plus and negative gain span */
  idxtype maxgain;
  ListNodeType *nodes;
  ListNodeType **buckets;

  /* Heap version of the data structure */
  KeyValueType *heap;
  idxtype *locator;
};

typedef struct PQueueType PQueueType;


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edegreedef {
  idxtype pid;
  idxtype ed;
};
typedef struct edegreedef EDegreeType;


/*************************************************************************
* The following data structure stores an edge for vol
**************************************************************************/
struct vedegreedef {
  idxtype pid;
  idxtype ed, ned;
  idxtype gv;
};
typedef struct vedegreedef VEDegreeType;


/*************************************************************************
* This data structure holds various working space data
**************************************************************************/
struct workspacedef {
  idxtype *core;			/* Where pairs, indices, and degrees are coming from */
  idxtype maxcore, ccore;

  EDegreeType *edegrees;
  VEDegreeType *vedegrees;
  idxtype cdegree;

  idxtype *auxcore;			/* This points to the memory of the edegrees */

  idxtype *pmat;			/* An array of k^2 used for eliminating domain 
                                           connectivity in k-way refinement */
};

typedef struct workspacedef WorkSpaceType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 idxtype id, ed;            	/* ID/ED of nodes */
 idxtype ndegrees;          	/* The number of different ext-degrees */
 EDegreeType *edegrees;     	/* List of edges */
};

typedef struct rinfodef RInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* vol-based partition
**************************************************************************/
struct vrinfodef {
 idxtype id, ed, nid;            	/* ID/ED of nodes */
 idxtype gv;            		/* IV/EV of nodes */
 idxtype ndegrees;          	/* The number of different ext-degrees */
 VEDegreeType *edegrees;     	/* List of edges */
};

typedef struct vrinfodef VRInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 idxtype edegrees[2];  
};

typedef struct nrinfodef NRInfoType;


/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  idxtype nvtxs, nedges;	/* The # of vertices and edges in the graph */
  idxtype *xadj;		/* Pointers to the locally stored vertices */
  idxtype *vwgt;		/* Vertex weights */
  idxtype *vsize;		/* Vertex sizes for min-volume formulation */
  idxtype *adjncy;		/* Array that stores the adjacency lists of nvtxs */
  idxtype *adjwgt;		/* Array that stores the weights of the adjacency lists */


  /* These are to keep track control if the corresponding fields correspond to
     application or library memory */
  int free_xadj, free_vwgt, free_vsize, free_adjncy, free_adjwgt;


  double *coords;               /* x, y, and z, Coordinates */

  idxtype *adjwgtsum;		/* The sum of the adjacency weight of each vertex */

  idxtype *label;

  idxtype *cmap;

  /* Partition parameters */
  idxtype mincut, minvol;
  idxtype *where, *pwgts;
  idxtype nbnd;
  idxtype *bndptr, *bndind;

  /* Bisection refinement parameters */
  idxtype *id, *ed;

  /* K-way refinement parameters */
  RInfoType *rinfo;

  /* K-way volume refinement parameters */
  VRInfoType *vrinfo;

  /* Node refinement information */
  NRInfoType *nrinfo;


  /* Additional info needed by the MOC routines */
  idxtype ncon;			/* The # of constrains */ 
  float *nvwgt;			/* Normalized vertex weights */
  float *npwgts;		/* The normalized partition weights */

  struct graphdef *coarser, *finer;
};

typedef struct graphdef GraphType;




/*************************************************************************
* The following structure stores information used by Metis
**************************************************************************/
struct controldef {
  idxtype CoarsenTo;		/* The # of vertices in the coarsest graph */
  idxtype dbglvl;			/* Controls the debuging output of the program */
  idxtype CType;			/* The type of coarsening */
  idxtype IType;			/* The type of initial partitioning */
  idxtype RType;			/* The type of refinement */
  idxtype maxvwgt;			/* The maximum allowed weight for a vertex */
  float nmaxvwgt;		/* The maximum allowed weight for a vertex for each constrain */
  idxtype optype;			/* Type of operation */
  idxtype pfactor;			/* .1*prunning factor */
  idxtype nseps;			/* The number of separators to be found during multiple bisections */
  idxtype oflags;

  WorkSpaceType wspace;		/* Work Space Informations */

  /* Various Timers */
  double TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, UncoarsenTmr, 
         SepTmr, RefTmr, ProjectTmr, SplitTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;

};

typedef struct controldef CtrlType;


/*************************************************************************
* The following data structure stores max-partition weight info for 
* Vertical MOC k-way refinement
**************************************************************************/
struct vpwgtdef {
  float max[2][MAXNCON];
  idxtype imax[2][MAXNCON];
};

typedef struct vpwgtdef VPInfoType;




