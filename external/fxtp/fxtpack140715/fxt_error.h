#ifndef __FXT_ERROR_H_INCLUDED__
#define __FXT_ERROR_H_INCLUDED__

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  error handling for fxtpack
*********************************************************************/

/*** error levels ***/

#define FXT_ERROR_CLEAR  0	/* no error, normal state */
#define FXT_ERROR_WARN   1	/* warning */
#define FXT_ERROR_USAGE  2	/* erroneous usage */
#define FXT_ERROR_FXTBUG 3	/* bug of the fxtpack */
#define FXT_ERROR_SYSTEM 4  	/* system error (malloc fail etc) */

#define FXT_ERROR_MSG_SIZE 1024 /* size of error message buffer */

/*** clear errors ***/
void fxt_error_clear(void);

/*** raise error level from USAGE to FXTBUG ***/
int fxt_error_raise(void);

/*** returns the error level ***/
int fxt_error_level(void);

/*** returns the error message: do not destruct returned string ***/
char* fxt_error_message(void);

/*** set error level and message ***/
void fxt_error_set(int level, char *fmt, ...);

#endif
