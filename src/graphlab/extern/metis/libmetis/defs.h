/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h,v 1.4 2003/01/29 15:54:41 karypis Exp $
 *
 */

#define METISTITLE              "  METIS 5.0 Copyright 2001-06, Regents of the University of Minnesota\n\n"
#define MAXLINE			1280000

#define LTERM			(void **) 0	/* List terminator for gk_free() */

#define MAXNCON			16		/* The maximum number of constrains */
#define MAXNOBJ			16		/* The maximum number of objectives */

#define PLUS_GAINSPAN   	500             /* Parameters for FM buckets */
#define NEG_GAINSPAN    	500

#define HTLENGTH		((1<<11)-1)

/* Meaning of various options[] parameters */
#define OPTION_PTYPE		0
#define OPTION_CTYPE		1
#define OPTION_ITYPE		2
#define OPTION_RTYPE		3
#define OPTION_DBGLVL		4
#define OPTION_OFLAGS		5
#define OPTION_PFACTOR		6
#define OPTION_NSEPS		7

#define OFLAG_COMPRESS		1	/* Try to compress the graph */
#define OFLAG_CCMP		2	/* Find and order connected components */


/* Default options for PMETIS */
#define PMETIS_CTYPE		MTYPE_SHEM
#define PMETIS_ITYPE		ITYPE_GGPKL
#define PMETIS_RTYPE		RTYPE_FM
#define PMETIS_DBGLVL		0

/* Default options for KMETIS */
#define KMETIS_CTYPE		MTYPE_SHEM
#define KMETIS_ITYPE		ITYPE_PMETIS
#define KMETIS_RTYPE		RTYPE_KWAYRANDOM_MCONN
#define KMETIS_DBGLVL		0

/* Default options for OEMETIS */
#define OEMETIS_CTYPE		MTYPE_SHEM
#define OEMETIS_ITYPE		ITYPE_GGPKL
#define OEMETIS_RTYPE		RTYPE_FM
#define OEMETIS_DBGLVL		0

/* Default options for ONMETIS */
#define ONMETIS_CTYPE		MTYPE_SHEM
#define ONMETIS_ITYPE		ITYPE_GGPKL
#define ONMETIS_RTYPE		RTYPE_SEP1SIDED
#define ONMETIS_DBGLVL		0
#define ONMETIS_OFLAGS		OFLAG_COMPRESS
#define ONMETIS_PFACTOR		-1
#define ONMETIS_NSEPS		1

/* Default options for McPMETIS */
#define McPMETIS_CTYPE		MTYPE_SHEBM_ONENORM
#define McPMETIS_ITYPE		ITYPE_RANDOM
#define McPMETIS_RTYPE		RTYPE_FM
#define McPMETIS_DBGLVL		0

/* Default options for McKMETIS */
#define McKMETIS_CTYPE		MTYPE_SHEBM_ONENORM
#define McKMETIS_ITYPE		ITYPE_McHPMETIS
#define McKMETIS_RTYPE		RTYPE_KWAYRANDOM
#define McKMETIS_DBGLVL		0

/* Default options for KVMETIS */
#define KVMETIS_CTYPE		MTYPE_SHEM
#define KVMETIS_ITYPE		ITYPE_PMETIS
#define KVMETIS_RTYPE		RTYPE_KWAYRANDOM
#define KVMETIS_DBGLVL		0


/* Operations supported by stand-alone code */
#define OP_PMETIS		1
#define OP_KMETIS		2
#define OP_OEMETIS		3
#define OP_ONMETIS		4
#define OP_ONWMETIS		5
#define OP_KVMETIS		6


#define UNMATCHED		-1

#define HTABLE_EMPTY    	-1

#define NGR_PASSES		4	/* Number of greedy refinement passes */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define LARGENIPARTS		8	/* Number of random initial partitions */
#define SMALLNIPARTS		3	/* Number of random initial partitions */

#define COARSEN_FRACTION	0.75	/* Node reduction between succesive coarsening levels */
#define COARSEN_FRACTION2	0.99	/* Node reduction between succesive coarsening levels */
#define UNBALANCE_FRACTION		1.03

#define COMPRESSION_FRACTION		0.85

#define ORDER_UNBALANCE_FRACTION	1.10

#define MMDSWITCH		        200

#define HORIZONTAL_IMBALANCE		1.05

/* Debug Levels */
#define DBG_TIME	1		/* Perform timing analysis */
#define DBG_OUTPUT	2
#define DBG_COARSEN   	4		/* Show the coarsening progress */
#define DBG_REFINE	8		/* Show info on communication during folding */
#define DBG_IPART	16		/* Show info on initial partition */
#define DBG_MOVEINFO	32		/* Show info on communication during folding */
#define DBG_KWAYPINFO	64		/* Show info on communication during folding */
#define DBG_SEPINFO	128		/* Show info on communication during folding */
