#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_sarld.h"
#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  creation and deletion of simple array representation
*********************************************************************/

fxt_sarld* fxt_sarld_new(long nrow, long ncol, long ntmp, long nent) {
  fxt_sarld *ar;

  ar = (fxt_sarld*) malloc(sizeof(fxt_sarld));
  if (ar == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_sarld_new: allocation failed\n");
    return NULL;
  }

  ar->nrow = nrow;
  ar->ncol = ncol;
  ar->ntmp = ntmp;

  ar->p = (long*) malloc(sizeof(long) * (nent + ntmp + ncol + 1));
  if (ar->p == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_sarld_new: allocation failed\n");
    return NULL;
  }

  ar->v = (double*) malloc(sizeof(double) * nent);
  if (ar->v == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_sarld_new: allocation failed\n");
    return NULL;
  }

  ar->k = ar->p + ntmp + ncol + 1;
  ar->p[0] = 0;

  return ar;
}


void fxt_sarld_del(fxt_sarld *ar) {

  /* check null pointers */
  if (ar == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_del: null pointer\n");
    return;
  }

  free(ar->v);
  ar->v = NULL;

  free(ar->p);
  ar->p = NULL;

  free(ar);
}


void fxt_sarld_save(FILE *f, fxt_sarld *ar) {
  long size[4];

  /* check null pointers */
  if (f == NULL || ar == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_save: null pointer\n");
    return;
  }

  size[0] = ar->nrow;
  size[1] = ar->ncol;
  size[2] = ar->ntmp;
  size[3] = ar->p[ar->ncol + ar->ntmp];

  fxt_file_writelongs(f, size, 4);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_file_writelongs(f, ar->p, ar->ncol + ar->ntmp + 1 + size[3]);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_file_writedoubles(f, ar->v, size[3]);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}


fxt_sarld* fxt_sarld_restore(FILE *f) {
  long size[4];
  fxt_sarld *ar;

  /* check null pointers */
  if (f == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_restore: null pointer\n");
    return NULL;
  }

  fxt_file_readlongs(f, size, 4);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  ar = fxt_sarld_new(size[0], size[1], size[2], size[3]);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_file_readlongs(f, ar->p, ar->ncol + ar->ntmp + 1 + size[3]);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_file_readdoubles(f, ar->v, size[3]);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  return ar;
}
