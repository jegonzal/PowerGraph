/*!
\file gk_struct.h
\brief This file contains various datastructures used/provided by GKlib

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_struct.h 1421 2007-04-06 14:37:41Z karypis $ \endverbatim
*/

#ifndef _GK_STRUCT_H_
#define _GK_STRUCT_H_


/********************************************************************/
/*! Generator for gk_??KeyVal_t data structure */
/********************************************************************/
#define GK_KEYVALUE_T(NAME, KEYTYPE, VALTYPE) \
typedef struct {\
  KEYTYPE key;\
  VALTYPE val;\
} NAME\

/* The actual KeyVal data structures */
GK_KEYVALUE_T(gk_ckv_t,   char,     gk_idx_t);
GK_KEYVALUE_T(gk_ikv_t,   int,      gk_idx_t);
GK_KEYVALUE_T(gk_fkv_t,   float,    gk_idx_t);
GK_KEYVALUE_T(gk_dkv_t,   double,   gk_idx_t);
GK_KEYVALUE_T(gk_skv_t,   char *,   gk_idx_t);
GK_KEYVALUE_T(gk_idxkv_t, gk_idx_t, gk_idx_t);





/*-------------------------------------------------------------
 * The following data structure stores a sparse CSR format
 *-------------------------------------------------------------*/
typedef struct {
  int nrows, ncols;
  int *rowptr, *colptr;
  int *rowind, *colind;
  int *irowval, *icolval;
  float *frowval, *fcolval;
} gk_CSR_t;


/*-------------------------------------------------------------
 * The following data structure implements max-heap priority
 * queue
 *-------------------------------------------------------------*/
typedef struct {
  int nnodes;
  int maxnodes;

  /* Heap version of the data structure */
  gk_fkv_t *heap;
  int *locator;
} gk_PQueue_t;




/*-------------------------------------------------------------
* The following data structure implements a string-2-int mapping
* table used for parsing command-line options
*-------------------------------------------------------------*/
typedef struct {
  char *name;
  int id;
} gk_StringMap_t;


/*------------------------------------------------------------
 * This structure implements a simple hash table
 *------------------------------------------------------------*/
typedef struct {
  int nelements;                /* The overall size of the hash-table */
  int htsize;                   /* The current size of the hash-table */
  gk_ikv_t *harray;       /* The actual hash-table */
} gk_HTable_t;


/*------------------------------------------------------------
 * This structure implements a gk_Tokens_t list returned by the
 * string tokenizer
 *------------------------------------------------------------*/
typedef struct {
  int ntoks;        /* The number of tokens in the input string */
  char *strbuf;     /* The memory that stores all the entries */
  char **list;      /* Pointers to the strbuf for each element */
} gk_Tokens_t;

/*------------------------------------------------------------
 * This structure implements storage for an atom in a pdb file
 *------------------------------------------------------------*/
typedef struct {
	int       serial;
	char      *name;
	char			altLoc;
	char      *resname;
	char			chainid;	
	int       rserial;
	char			icode;
  char      element;
	double    x;
	double    y;
	double    z;
	double    opcy;
	double    tmpt;
} atom;

/*------------------------------------------------------------
 * This structure implements storage for a pdb protein 
 *------------------------------------------------------------*/
typedef struct {
	int natoms;			/* Number of atoms */
	int nresidues;  /* Number of residues based on coordinates */
	int ncas;
	int nbbs;
	int corruption;
	char *resSeq;	  /* Residue sequence based on coordinates    */
	atom *atoms;
	atom **bbs;
	atom **cas;
} pdbf;



/*************************************************************
 * Localization Structures for converting characters to integers
 *************************************************************/
typedef struct {
    int n;
    char *i2c;
    int *c2i;
} gk_i2cc2i_t;
 

/*******************************************************************
 *This structure implements storage of a protein sequence
 * *****************************************************************/
typedef struct {
    
    int len; /*Number of Residues */
    int *sequence; /* Stores the sequence*/
    
    
    int **pssm; /* Stores the pssm matrix */
    int **psfm; /* Stores the psfm matrix */
    char *name; /* Stores the name of the sequence */

    int nsymbols;

    
} gk_seq_t;


#endif
