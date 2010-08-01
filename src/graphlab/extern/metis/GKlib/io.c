/*!
\file  io.c
\brief Various file I/O functions.

This file contains various functions that perform I/O.

\date Started 4/10/95
\author George
\version\verbatim $Id: io.c 1430 2007-04-07 17:53:07Z karypis $ \endverbatim
*/


#include <GKlib.h>

/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *gk_fopen(char *fname, char *mode, const char *msg)
{
  FILE *fp;
  char errmsg[8192];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(errmsg);
  errexit("Failed on gk_fopen()\n");

  return NULL;
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void gk_fclose(FILE *fp)
{
  fclose(fp);
}



/*************************************************************************
* This function is the GKlib implementation of glibc's getline()
* function.
**************************************************************************/
gk_loop_t gk_getline(char **lineptr, size_t *n, FILE *stream)
{
  size_t i;
  int ch;

  if (feof(stream))
    return -1;  

  /* Initial memory allocation if *lineptr is NULL */
  if (*lineptr == NULL || *n == 0) {
    *n = 1024;
    *lineptr = gk_malloc((*n)*sizeof(char), "gk_getline: lineptr");
  }

  /* get into the main loop */
  for (i = 0; (ch = getc(stream)) != EOF; i++) {
    (*lineptr)[i] = (char)ch;

    /* reallocate memory if reached at the end of the buffer. The +2 is for '\0' */
    if (i+2 == *n) { 
      *n = 2* (*n);
      *lineptr = gk_realloc(*lineptr, (*n)*sizeof(char), "gk_getline: lineptr");
    }
      
    if (ch == '\n')
      break;
  }
  (*lineptr)[i+1] = '\0';

  return i;
}

