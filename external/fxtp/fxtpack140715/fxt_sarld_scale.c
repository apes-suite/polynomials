#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_sarld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  scaling for simple array representation
*********************************************************************/

void fxt_sarld_scale(fxt_sarld *ar, double s) {
  long j;

  /* check null pointers */
  if (ar == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_scale: null pointer\n");
    return;
  }

  /* scale */
  for (j=0; j< ar->p[ar->ncol]; j++)
    ar->v[j] *= s;
}

void fxt_sarld_scalerow(fxt_sarld *ar, fxt_vecld *s) {
  long j;

  /* check null pointers */
  if (ar == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_scalerow: null pointer\n");
    return;
  }

  /* check size */
  if (s->n != ar->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_scalerow: size mismatch\n");
    return;
  }

  /* scale rows */
  for (j=0; j< ar->p[ar->ncol + ar->ntmp]; j++) {
    long k = ar->k[j] - ar->ncol;
    if (0 <= k && k < ar->nrow)
      ar->v[j] *= s->v[k];
  }
}
