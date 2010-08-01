/*!
\file gk_blas.h
\brief This file contains GKlib's code generators for BLAS-like routines

\date   Started 3/28/07
\author George
\version\verbatim $Id: gk_blas.h 1421 2007-04-06 14:37:41Z karypis $ \endverbatim
*/

#ifndef _GK_BLAS_H_
#define _GK_BLAS_H_


/*************************************************************************/
/*! The macro for gk_?set()-class of routines */
/*************************************************************************/
#define GK_SET(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE val, TYPE *x)\
{\
  size_t i;\
\
  for (i=0; i<n; i++)\
    x[i] = val;\
\
  return x;\
}

#define GK_SET_PROTO(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE val, TYPE *x);




/*************************************************************************/
/*! The macro for gk_?incset()-class of routines */
/*************************************************************************/
#define GK_INCSET(NAME, TYPE)\
TYPE *NAME(size_t n, TYPE baseval, TYPE *x)\
{\
  size_t i;\
\
  for (i=0; i<n; i++)\
    x[i] = baseval+i;\
\
  return x;\
}

#define GK_INCSET_PROTO(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE baseval, TYPE *x);



/*************************************************************************/
/*! The macro for gk_?argmax()-class of routines */
/*************************************************************************/
#define GK_ARGMAX(NAME, TYPE) \
size_t NAME(size_t n, TYPE *x)\
{\
  size_t i, max=0;\
\
  for (i=1; i<n; i++)\
    max = (x[i] > x[max] ? i : max);\
\
  return max;\
}

#define GK_ARGMAX_PROTO(NAME, TYPE) \
size_t NAME(size_t n, TYPE *x);



/*************************************************************************/
/*! The macro for gk_?argmin()-class of routines */
/*************************************************************************/
#define GK_ARGMIN(NAME, TYPE) \
size_t NAME(size_t n, TYPE *x)\
{\
  size_t i, min=0;\
\
  for (i=1; i<n; i++)\
    min = (x[i] < x[min] ? i : min);\
\
  return min;\
}

#define GK_ARGMIN_PROTO(NAME, TYPE) \
size_t NAME(size_t n, TYPE *x);





/*************************************************************************/
/*! The macro for gk_?argmax_n()-class of routines */
/*************************************************************************/
#define GK_ARGMAX_N(NAME, TYPE, KEYVAL_T, KEYVALMALLOC, KEYVALSORT) \
size_t NAME(size_t n, TYPE *x, size_t k)\
{\
  size_t i, max_n;\
  KEYVAL_T *cand;\
\
  cand = KEYVALMALLOC(n, "GK_ARGMAX_N: cand");\
\
  for (i=0; i<n; i++) {\
    cand[i].val = i;\
    cand[i].key = x[i];\
  }\
  KEYVALSORT(n, cand);\
\
  max_n = cand[k-1].val;\
\
  gk_free((void *)&cand, LTERM);\
\
  return max_n;\
}

#define GK_ARGMAX_N_PROTO(NAME, TYPE) \
size_t NAME(size_t n, TYPE *x, size_t k);




/*************************************************************************/
/*! The macro for gk_?sum()-class of routines */
/**************************************************************************/
#define GK_SUM(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx)\
{\
  size_t i;\
  OUTTYPE sum = 0;\
\
  for (i=0; i<n; i++, x+=incx)\
    sum += (*x);\
\
  return sum;\
}

#define GK_SUM_PROTO(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx);





/*************************************************************************/
/*! The macro for gk_?scale()-class of routines */
/**************************************************************************/
#define GK_SCALE(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE alpha, TYPE *x, size_t incx)\
{\
  size_t i;\
\
  for (i=0; i<n; i++, x+=incx)\
    (*x) *= alpha;\
\
  return x;\
}

#define GK_SCALE_PROTO(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE alpha, TYPE *x, size_t incx);




/*************************************************************************/
/*! The macro for gk_?norm2()-class of routines */
/**************************************************************************/
#define GK_NORM2(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx)\
{\
  size_t i;\
  OUTTYPE partial = 0;\
\
  for (i=0; i<n; i++, x+=incx)\
    partial += (*x) * (*x);\
\
  return (partial > 0 ? (OUTTYPE)sqrt((double)partial) : 0.0);\
}

#define GK_NORM2_PROTO(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx);




/*************************************************************************/
/*! The macro for gk_?dot()-class of routines */
/**************************************************************************/
#define GK_DOT(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx, INTYPE *y, size_t incy)\
{\
  size_t i;\
  OUTTYPE partial = 0.0;\
 \
  for (i=0; i<n; i++, x+=incx, y+=incy)\
    partial += (*x) * (*y);\
\
  return partial;\
}

#define GK_DOT_PROTO(NAME, INTYPE, OUTTYPE) \
OUTTYPE NAME(size_t n, INTYPE *x, size_t incx, INTYPE *y, size_t incy);



/*************************************************************************/
/*! The macro for gk_?axpy()-class of routines */
/**************************************************************************/
#define GK_AXPY(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE alpha, TYPE *x, size_t incx, TYPE *y, size_t incy)\
{\
  size_t i;\
  TYPE *y_in = y;\
\
  for (i=0; i<n; i++, x+=incx, y+=incy)\
    *y += alpha*(*x);\
\
  return y_in;\
}

#define GK_AXPY_PROTO(NAME, TYPE) \
TYPE *NAME(size_t n, TYPE alpha, TYPE *x, size_t incx, TYPE *y, size_t incy);

#endif
