/*!
\file  gk_memory.h
\brief This file contains GKlib's code generators for memory allocation routines

\date   Started 3/29/07
\author George
\version\verbatim $Id: gk_memory.h 1432 2007-04-07 20:06:19Z karypis $ \endverbatim
*/

#ifndef _GK_MEMORY_H_
#define _GK_MEMORY_H_


/*************************************************************************/
/*! The macro for gk_?malloc()-class of routines */
/**************************************************************************/
#define GK_XMALLOC(NAME, TYPE)\
TYPE *NAME(size_t n, char *msg)\
{ return (TYPE *)gk_malloc(sizeof(TYPE)*n, msg); }

#define GK_XMALLOC_PROTO(NAME, TYPE)\
TYPE *NAME(size_t n, char *msg);


/*************************************************************************/
/*! The macro for gk_?realloc()-class of routines */
/**************************************************************************/
#define GK_XREALLOC(NAME, TYPE)\
TYPE *NAME(TYPE *ptr, size_t n, char *msg)\
{ return (TYPE *)gk_realloc((void *)ptr, sizeof(TYPE)*n, msg); }

#define GK_XREALLOC_PROTO(NAME, TYPE)\
TYPE *NAME(TYPE *ptr, size_t n, char *msg);


/*************************************************************************/
/*! The macro for gk_?smalloc()-class of routines */
/**************************************************************************/
#define GK_XSMALLOC(NAME, TYPE, SET)\
TYPE *NAME(size_t n, TYPE ival, char *msg)\
{ return SET(n, ival, (TYPE *)gk_malloc(sizeof(TYPE)*n, msg)); }

#define GK_XSMALLOC_PROTO(NAME, TYPE)\
TYPE *NAME(size_t n, TYPE ival, char *msg);




/*************************************************************************/
/*! The macro for gk_?AllocMatrix()-class of routines */
/**************************************************************************/
#define GK_ALLOCMATRIX(NAME, TYPE, ROWALLOC)\
TYPE **NAME(size_t ndim1, size_t ndim2, TYPE value, char *errmsg)\
{\
  gk_loop_t i;\
  TYPE **matrix;\
\
  matrix = (TYPE **)gk_malloc(ndim1*sizeof(TYPE *), errmsg);\
  for (i=0; i<ndim1; i++) \
    matrix[i] = ROWALLOC(ndim2, value, errmsg);\
\
  return matrix;\
}

#define GK_ALLOCMATRIX_PROTO(NAME, TYPE) \
TYPE **NAME(size_t ndim1, size_t ndim2, TYPE value, char *errmsg);


/*************************************************************************/
/*! The macro for gk_?AllocMatrix()-class of routines */
/**************************************************************************/
#define GK_FREEMATRIX(NAME, TYPE)\
void NAME(TYPE ***r_matrix, size_t ndim1, size_t ndim2)\
{\
  gk_loop_t i;\
  TYPE **matrix;\
\
  matrix = *r_matrix;\
\
  for (i=0; i<ndim1; i++) \
    gk_free((void **)&(matrix[i]), LTERM);\
\
  gk_free((void **)matrix, LTERM);\
\
  *r_matrix = NULL;\
}

#define GK_FREEMATRIX_PROTO(NAME, TYPE)\
void NAME(TYPE ***r_matrix, size_t ndim1, size_t ndim2);



/*************************************************************************/
/*! The macro for gk_?SetMatrix()-class of routines */
/**************************************************************************/
#define GK_SETMATRIX(NAME, TYPE)\
void NAME(TYPE **matrix, size_t ndim1, size_t ndim2, TYPE value)\
{\
  gk_loop_t i, j;\
\
  for (i=0; i<ndim1; i++) {\
    for (j=0; j<ndim2; j++)\
      matrix[i][j] = value;\
  }\
}

#define GK_SETMATRIX_PROTO(NAME, TYPE)\
void NAME(TYPE **matrix, size_t ndim1, size_t ndim2, TYPE value);

#endif
