/*!
\file gk_defs.h
\brief This file contains various constants definitions

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_defs.h 1277 2007-03-27 21:17:33Z karypis $ \endverbatim
*/

#ifndef _GK_DEFS_H_
#define _GK_DEFS_H_


#define LTERM                   (void **) 0     /* List terminator for GKfree() */

#define HTABLE_EMPTY            -1
#define HTABLE_DELETED          -2
#define HTABLE_FIRST             1
#define HTABLE_NEXT              2

/* pdb corruption bit switches */
#define CRP_ALTLOCS    1
#define CRP_MISSINGCA  2
#define CRP_MISSINGBB  4
#define CRP_MULTICHAIN 8
#define CRP_MULTICA    16
#define CRP_MULTIBB    32

#define MAXLINELEN 300000

/* custom signals */
#define SIGMEM                  SIGABRT
#define SIGERR                  SIGABRT

#endif
