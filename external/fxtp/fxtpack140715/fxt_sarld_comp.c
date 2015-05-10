#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_sarld.h"

#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  computations of simple array representation
*********************************************************************/

void fxt_sarld_evl(fxt_sarld *ar, double *w) {
  long i, j;  double t;

  /* check null pointers */
  if (ar == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_evl: null pointer\n");
    return;
  }

  for (i= ar->ncol; i < ar->ncol + ar->nrow + ar->ntmp; i++)
    w[i] = 0.0;

  for (i= ar->ncol - 1; i >= 0; i--) {
    t = w[i];
    for (j= ar->p[i]; j< ar->p[i+1]; j++)
      w[ar->k[j]] += ar->v[j] * t;
  }

  for (i= ar->ncol + ar->ntmp - 1; i >= ar->ncol; i--) {
    t = w[ar->nrow + i];
    for (j= ar->p[i]; j< ar->p[i+1]; j++)
      w[ar->k[j]] += ar->v[j] * t;
  }
}

void fxt_sarld_exp(fxt_sarld *ar, double *w) {
  long i, j;  double t;

  /* check null pointers */
  if (ar == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_exp: null pointer\n");
    return;
  }

  for (i= ar->ncol; i< ar->ncol + ar->ntmp; i++) {
    t = 0.0;
    for (j= ar->p[i]; j< ar->p[i+1]; j++)
      t += ar->v[j] * w[ar->k[j]];
    w[ar->nrow + i] = t;
  }

  for (i= 0; i< ar->ncol; i++) {
    t = 0.0;
    for (j= ar->p[i]; j< ar->p[i+1]; j++)
      t += ar->v[j] * w[ar->k[j]];
    w[i] = t;
  }
}
