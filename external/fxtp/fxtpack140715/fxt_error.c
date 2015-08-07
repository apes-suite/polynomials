#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "fxt_error.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  error handling for FXTPACK
*********************************************************************/


static int error_level = 0;
static char error_message[FXT_ERROR_MSG_SIZE];


/*** clear errors ***/
void fxt_error_clear(void) {
  error_level = 0;
  error_message[0] = '\0';
}


/*** raise error level from USAGE to FXTBUG ***/
int fxt_error_raise(void) {
  if (error_level == FXT_ERROR_USAGE)
    error_level = FXT_ERROR_FXTBUG;

  return error_level;
}


/*** returns the error level ***/
int fxt_error_level(void) {
  return error_level;
}


/*** returns the error message: do not distruct it! ***/
char* fxt_error_message(void) {
  return error_message;
}


/*** set error level and message: %d, %f, %e, %s available ***/
void fxt_error_set(int level, char *fmt, ...) {
  va_list ap;			/* arguments */
  char *p, *s;

  /* ignore lower level errors */
  if (level < error_level)
    return;

  /* set new error level */
  error_level = level;

  /* initialize ap */
  va_start(ap, fmt);

  /* initialize write pointer */
  s = &(error_message[0]);

  /* read format */
  for (p = fmt; *p != '\0'; p++) {
    if (*p != '%')		/* usual character */
      *(s++) = *p;
    else {
      switch (p[1]) {		/* in case of '%' */
      case 'd':
	{			/* integer */
	  int val;
	  val = va_arg(ap, int);
	  sprintf(s, "%d", val);
	}
	break;

      case 'f':  case 'e':
	{			/* double */
	  double val;
	  val = va_arg(ap, double);
	  sprintf(s, (p[1] == 'f' ? "%f" : "%e"), val);
	}
	break;

      case 's':
	{			/* string */
	  char *val;
	  val = va_arg(ap, char*);
	  sprintf(s, "%s", val);
	}
	break;

      default:			/* just ignore '%' */
	continue;
      }

      /* skip d, f, e or s */
      p ++;

      /* advance write pointer */
      while (*s != '\0')
	s ++;
    }
  }

  /* terminate message string */
  *s = '\0';

  /* finalize argument list */
  va_end(ap);
}
