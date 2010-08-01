/*!
\file blas.c
\brief This file contains GKlib's implementation of BLAS-like routines

The BLAS routines that are currently implemented are mostly level-one.
They follow a naming convention of the type gk_[type][name], where
[type] is one of c, i, f, and d, based on C's four standard scalar
datatypes of characters, integers, floats, and doubles.

These routines are implemented using a generic macro template,
which is used for code generation.

\date   Started 9/28/95
\author George
\version\verbatim $Id: blas.c 1421 2007-04-06 14:37:41Z karypis $ \endverbatim
*/

#include <GKlib.h>



/*************************************************************************/
/*! Set the elements of an array to a specified value */
/*************************************************************************/
GK_SET(gk_cset,   char)
GK_SET(gk_iset,   int)
GK_SET(gk_fset,   float)
GK_SET(gk_dset,   double)
GK_SET(gk_idxset, gk_idx_t)



/*************************************************************************/
/*! Set the elements of an array in increasing order */
/*************************************************************************/
GK_INCSET(gk_cincset,   char)
GK_INCSET(gk_iincset,   int)
GK_INCSET(gk_fincset,   float)
GK_INCSET(gk_dincset,   double)
GK_INCSET(gk_idxincset, gk_idx_t)



/*************************************************************************/
/*! Return the index of the maximum element in a vector */
/*************************************************************************/
GK_ARGMAX(gk_cargmax,   char)
GK_ARGMAX(gk_iargmax,   int)
GK_ARGMAX(gk_fargmax,   float)
GK_ARGMAX(gk_dargmax,   double)
GK_ARGMAX(gk_idxargmax, gk_idx_t)



/*************************************************************************/
/*! Return the index of the minimum element in a vector */
/*************************************************************************/
GK_ARGMIN(gk_cargmin,   char)
GK_ARGMIN(gk_iargmin,   int)
GK_ARGMIN(gk_fargmin,   float)
GK_ARGMIN(gk_dargmin,   double)
GK_ARGMIN(gk_idxargmin, gk_idx_t)



/*************************************************************************/
/*! Return the index of the Kth largest element in a vector */
/*************************************************************************/
GK_ARGMAX_N(gk_cargmax_n,   char,     gk_ckv_t,   gk_ckvmalloc,   gk_dckvsort)
GK_ARGMAX_N(gk_iargmax_n,   int,      gk_ikv_t,   gk_ikvmalloc,   gk_dikvsort)
GK_ARGMAX_N(gk_fargmax_n,   float,    gk_fkv_t,   gk_fkvmalloc,   gk_dfkvsort)
GK_ARGMAX_N(gk_dargmax_n,   double,   gk_dkv_t,   gk_dkvmalloc,   gk_ddkvsort)
GK_ARGMAX_N(gk_idxargmax_n, gk_idx_t, gk_idxkv_t, gk_idxkvmalloc, gk_didxkvsort)


/*************************************************************************/
/*! Sums the entries of an array */
/*************************************************************************/
GK_SUM(gk_csum,   char,     intmax_t)
GK_SUM(gk_isum,   int,      intmax_t)
GK_SUM(gk_fsum,   float,    float)
GK_SUM(gk_dsum,   double,   double)
GK_SUM(gk_idxsum, gk_idx_t, intmax_t)



/*************************************************************************/
/*! Scales the elements of a vector by a constant */
/*************************************************************************/
GK_SCALE(gk_cscale,   char)
GK_SCALE(gk_iscale,   int)
GK_SCALE(gk_fscale,   float)
GK_SCALE(gk_dscale,   double)
GK_SCALE(gk_idxscale, gk_idx_t)


/*************************************************************************/
/*! Computes the 2-norm of a vector */
/*************************************************************************/
GK_NORM2(gk_cnorm2,   char,     intmax_t)
GK_NORM2(gk_inorm2,   int,      intmax_t)
GK_NORM2(gk_fnorm2,   float,    float)
GK_NORM2(gk_dnorm2,   double,   double)
GK_NORM2(gk_idxnorm2, gk_idx_t, intmax_t)




/*************************************************************************/
/*! Computes the dot-product of two vectors */
/*************************************************************************/
GK_DOT(gk_cdot,   char,     intmax_t)
GK_DOT(gk_idot,   int,      intmax_t)
GK_DOT(gk_fdot,   float,    float)
GK_DOT(gk_ddot,   double,   double)
GK_DOT(gk_idxdot, gk_idx_t, intmax_t)



/*************************************************************************/
/*! Computes y = ax+y */
/*************************************************************************/
GK_AXPY(gk_caxpy,   char)
GK_AXPY(gk_iaxpy,   int)
GK_AXPY(gk_faxpy,   float)
GK_AXPY(gk_daxpy,   double)
GK_AXPY(gk_idxaxpy, gk_idx_t)


