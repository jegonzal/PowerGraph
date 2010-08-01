/*!
\file gk_externs.h
\brief This file contains definitions of external variables created by GKlib

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_externs.h 1277 2007-03-27 21:17:33Z karypis $ \endverbatim
*/

#ifndef _GK_EXTERNS_H_
#define _GK_EXTERNS_H_

#include <metis.h>

/*************************************************************************
* Extern variable definition. Hopefully, the __thread makes them thread-safe.
**************************************************************************/
#ifndef _GK_ERROR_C_

#ifdef __linux__
extern __thread jmp_buf gk_return_to_entry; /* declared in error.c */
#else
extern jmp_buf gk_return_to_entry; /* declared in error.c */
#endif

#endif

#endif
