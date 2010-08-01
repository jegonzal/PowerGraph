/*!
\file  sort.c
\brief This file contains GKlib's various sorting routines

These routines are implemented using the GKSORT macro that is defined
in gk_qsort.h and is based on GNU's GLIBC qsort() implementation.

Additional sorting routines can be created using the same way that
these routines where defined.

\date   Started 4/4/07
\author George
\version\verbatim $Id: sort.c 1421 2007-04-06 14:37:41Z karypis $ \endverbatim
*/

#include <GKlib.h>



/*************************************************************************/
/*! Sorts an array of chars in increasing order */
/*************************************************************************/
void gk_icsort(int n, char *base)
{
#define char_lt(a, b) ((*a) < (*b))
  GKQSORT(char, base, n, char_lt);
}


/*************************************************************************/
/*! Sorts an array of chars in decreasing order */
/*************************************************************************/
void gk_dcsort(int n, char *base)
{
#define char_gt(a, b) ((*a) > (*b))
  GKQSORT(char, base, n, char_gt);
}


/*************************************************************************/
/*! Sorts an array of integers in increasing order */
/*************************************************************************/
void gk_iisort(int n, int *base)
{
#define int_lt(a, b) ((*a) < (*b))
  GKQSORT(int, base, n, int_lt);
}


/*************************************************************************/
/*! Sorts an array of integers in decreasing order */
/*************************************************************************/
void gk_disort(int n, int *base)
{
#define int_gt(a, b) ((*a) > (*b))
  GKQSORT(int, base, n, int_gt);
}


/*************************************************************************/
/*! Sorts an array of floats in increasing order */
/*************************************************************************/
void gk_ifsort(int n, float *base)
{
#define float_lt(a, b) ((*a) < (*b))
  GKQSORT(float, base, n, float_lt);
}


/*************************************************************************/
/*! Sorts an array of floats in decreasing order */
/*************************************************************************/
void gk_dfsort(int n, float *base)
{
#define float_gt(a, b) ((*a) > (*b))
  GKQSORT(float, base, n, float_gt);
}


/*************************************************************************/
/*! Sorts an array of doubles in increasing order */
/*************************************************************************/
void gk_idsort(int n, double *base)
{
#define double_lt(a, b) ((*a) < (*b))
  GKQSORT(double, base, n, double_lt);
}


/*************************************************************************/
/*! Sorts an array of doubles in decreasing order */
/*************************************************************************/
void gk_ddsort(int n, double *base)
{
#define double_gt(a, b) ((*a) > (*b))
  GKQSORT(double, base, n, double_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_idx_t in increasing order */
/*************************************************************************/
void gk_iidxsort(int n, gk_idx_t *base)
{
#define idx_lt(a, b) ((*a) < (*b))
  GKQSORT(gk_idx_t, base, n, idx_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_idx_t in decreasing order */
/*************************************************************************/
void gk_didxsort(int n, gk_idx_t *base)
{
#define idx_gt(a, b) ((*a) > (*b))
  GKQSORT(gk_idx_t, base, n, idx_gt);
}




/*************************************************************************/
/*! Sorts an array of gk_ckv_t in increasing order */
/*************************************************************************/
void gk_ickvsort(int n, gk_ckv_t *base)
{
#define ckey_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(gk_ckv_t, base, n, ckey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_ckv_t in decreasing order */
/*************************************************************************/
void gk_dckvsort(int n, gk_ckv_t *base)
{
#define ckey_gt(a, b) ((a)->key > (b)->key)
  GKQSORT(gk_ckv_t, base, n, ckey_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_ikv_t in increasing order */
/*************************************************************************/
void gk_iikvsort(int n, gk_ikv_t *base)
{
#define ikey_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(gk_ikv_t, base, n, ikey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_ikv_t in decreasing order */
/*************************************************************************/
void gk_dikvsort(int n, gk_ikv_t *base)
{
#define ikey_gt(a, b) ((a)->key > (b)->key)
  GKQSORT(gk_ikv_t, base, n, ikey_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_fkv_t in increasing order */
/*************************************************************************/
void gk_ifkvsort(int n, gk_fkv_t *base)
{
#define fkey_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(gk_fkv_t, base, n, fkey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_fkv_t in decreasing order */
/*************************************************************************/
void gk_dfkvsort(int n, gk_fkv_t *base)
{
#define fkey_gt(a, b) ((a)->key > (b)->key)
  GKQSORT(gk_fkv_t, base, n, fkey_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_dkv_t in increasing order */
/*************************************************************************/
void gk_idkvsort(int n, gk_dkv_t *base)
{
#define dkey_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(gk_dkv_t, base, n, dkey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_fkv_t in decreasing order */
/*************************************************************************/
void gk_ddkvsort(int n, gk_dkv_t *base)
{
#define dkey_gt(a, b) ((a)->key > (b)->key)
  GKQSORT(gk_dkv_t, base, n, dkey_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_skv_t in increasing order */
/*************************************************************************/
void gk_iskvsort(int n, gk_skv_t *base)
{
#define skey_lt(a, b) (strcmp((a)->key, (b)->key) < 0)
  GKQSORT(gk_skv_t, base, n, skey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_skv_t in decreasing order */
/*************************************************************************/
void gk_dskvsort(int n, gk_skv_t *base)
{
#define skey_gt(a, b) (strcmp((a)->key, (b)->key) > 0)
  GKQSORT(gk_skv_t, base, n, skey_gt);
}


/*************************************************************************/
/*! Sorts an array of gk_idxkv_t in increasing order */
/*************************************************************************/
void gk_iidxkvsort(int n, gk_idxkv_t *base)
{
#define idxkey_lt(a, b) ((a)->key < (b)->key)
  GKQSORT(gk_idxkv_t, base, n, idxkey_lt);
}


/*************************************************************************/
/*! Sorts an array of gk_idxkv_t in decreasing order */
/*************************************************************************/
void gk_didxkvsort(int n, gk_idxkv_t *base)
{
#define idxkey_gt(a, b) ((a)->key > (b)->key)
  GKQSORT(gk_idxkv_t, base, n, idxkey_gt);
}
