#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test error handling routines for FXTPACK
*********************************************************************/

void test_error(void) {
  char msg[FXT_ERROR_MSG_SIZE];

  printf("- test_error...\n");

  /* simple string */
  fxt_error_set(FXT_ERROR_WARN, "string");
  sprintf(msg, "string");

  if (fxt_error_level() != FXT_ERROR_WARN) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for warn\n");
    return;
  }

  if (strcmp(msg, fxt_error_message()) != 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: message mismatch\n");
    return;
  }

  /* check integer */
  fxt_error_set(FXT_ERROR_USAGE, "int %d val", 10);
  sprintf(msg, "int %d val", 10);

  if (fxt_error_level() != FXT_ERROR_USAGE) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for usage\n");
    return;
  }

  if (strcmp(msg, fxt_error_message()) != 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: message mismatch for %%d\n");
    return;
  }

  /* check raise */
  fxt_error_raise();
  if (fxt_error_level() != FXT_ERROR_FXTBUG) {
    fxt_error_clear();
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for raise\n");
    return;
  }

  /* double */
  fxt_error_set(FXT_ERROR_FXTBUG, "double %f %e", 0.5, 0.5);
  sprintf(msg, "double %f %e", 0.5, 0.5);

  if (fxt_error_level() != FXT_ERROR_FXTBUG) {
    fxt_error_clear();
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for bug\n");
    return;
  }

  if (strcmp(msg, fxt_error_message()) != 0) {
    fxt_error_clear();
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: message mismatch for %%f %%e\n");
    return;
  }

  /* check string */
  fxt_error_set(FXT_ERROR_SYSTEM, "str %s str", "string");
  sprintf(msg, "str %s str", "string");

  if (fxt_error_level() != FXT_ERROR_SYSTEM) {
    fxt_error_clear();
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for system\n");
    return;
  }

  if (strcmp(msg, fxt_error_message()) != 0) {
    fxt_error_clear();
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: message mismatch for %%s\n");
    return;
  }

  /* check clear */
  fxt_error_clear();
  
  if (fxt_error_level() != FXT_ERROR_CLEAR) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: level mismatch for clear\n");
    return;
  }

  if (strcmp("", fxt_error_message()) != 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_error: message mismatch for clear\n");
    return;
  }
}

/*
  returning errors of fxt_error_set via fxt_error_set has
  logical difficulties ...
*/
