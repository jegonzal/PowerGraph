/*
 * Copyright 2007, Regents of the University of Minnesota
 *
 * streamio.c
 *
 * This function contains various functions that deal with extending the
 * functionality of prinf and mscanf family of routines to understand
 * additional format conversions.
 *
 * At this point the following conversions are support:
 *  %D for idxtype (which either 32 or 64 bits)
 *
 * Started 4/6/07
 * George
 *
 * $Id: util.c,v 1.4 2003/04/13 04:45:12 karypis Exp $
 */

#include <metislib.h>


#define D_PATTERN               "((^%|[^%]%)[0-\\s\\+]*[0-9]*)D"
#define D_SCAN_REPLACEMENT      "$1" SCNIDX
#define D_PRINT_REPLACEMENT     "$1" PRIIDX



/*************************************************************************
* Custom mprintf() routine
**************************************************************************/
int mprintf(char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_PRINT_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vprintf((char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}


/*************************************************************************
* Custom msprintf() routine
**************************************************************************/
int msprintf(char *str, char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_PRINT_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vsprintf(str, (char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}

/*************************************************************************
* Custom mfprintf() routine
**************************************************************************/
int mfprintf(FILE *stream, char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_PRINT_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vfprintf(stream, (char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}


/*************************************************************************
* Custom mscanf() routine
**************************************************************************/
int mscanf(char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_SCAN_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vscanf((char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}


/*************************************************************************
* Custom msscanf() routine
**************************************************************************/
int msscanf(char *str, char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_SCAN_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vsscanf(str, (char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}

/*************************************************************************
* Custom mfscanf() routine
**************************************************************************/
int mfscanf(FILE *stream, char *format,...)
{
  char *new_format;
  int rc;
  va_list argp;

  gk_strstr_replace(format, D_PATTERN, D_SCAN_REPLACEMENT, "g", &new_format);

  /*mprintf("new_format: %s\n", new_format);*/

  va_start(argp, format);
  rc = vfscanf(stream, (char *)new_format, argp);
  va_end(argp);

  gk_free((void **)&new_format, LTERM);

  return rc;
}


