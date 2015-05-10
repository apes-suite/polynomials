#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_flptld.h"
#include "fxt_flptld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  Fast Legendre Polynomial Transform
*********************************************************************/

void fxt_flptld_evl(fxt_vecld *v, fxt_flptld *flpt,
		    fxt_vecld *u, fxt_vecld *w) {
  long npi, npx, i;
  fxt_sarld *ar;
  double t;

  /* check pointers */
  if (v == NULL || flpt == NULL || u == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_evl: null pointer\n");
    return;
  }

  /* check size */
  if (flpt->p != v->n || flpt->n + 1 < u->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_evl: size mismatch\n");
    return;
  }

  /* half-number of points, including zero */
  npi = (flpt->p + 1) / 2;

  /* half-number of points, excluding zero */
  npx = flpt->p / 2;

  /*** even transform ***/
  ar = flpt->ear;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_evl: insufficient working array\n");
    return;
  }

  /* copy even part of u to working array */
  for (i=0; i< (u->n + 1) / 2; i++)
    w->v[i] = u->v[2 * i];

  for (; i< ar->ncol; i++)
    w->v[i] = 0.0;

  /* do the even transform */
  fxt_sarld_evl(ar, w->v);

  /* store results to v from working array */

  /* treat zero specially */
  if (flpt->p % 2 == 1)
    v->v[npi] = w->v[ar->ncol + npi];

  /* copy other values */
  for (i=0; i< npx; i++) {
    t = w->v[ar->ncol + i];
    v->v[i] = t;
    v->v[flpt->p - i - 1] = t;
  }

  /*=== odd transform ===*/
  if (flpt->n == 0)
    return;

  ar = flpt->oar;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_evl: insufficient working array\n");
    return;
  }

  /* copy odd part of u to working array */
  for (i=0; i< u->n / 2; i++)
    w->v[i] = u->v[2 * i + 1];

  for (; i < ar->ncol; i++)
    w->v[i] = 0.0;

  /* do the odd transform */
  fxt_sarld_evl(ar, w->v);

  /* store results to v from working array */
  for (i=0; i< npx; i++) {
    t = w->v[ar->ncol + i];
    v->v[i] += t;
    v->v[flpt->p - i - 1] -= t;
  }
}

void fxt_flptld_exp(fxt_vecld *u, fxt_flptld *flpt,
		    fxt_vecld *v, fxt_vecld *w) {
  long npi, npx, i;
  fxt_sarld *ar;

  /* check pointers */
  if (v == NULL || flpt == NULL || u == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_exp: null pointer\n");
    return;
  }

  /* check size */
  if (flpt->p != v->n || flpt->n + 1 < u->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_exp: size mismatch\n");
    return;
  }

  /* half-number of points, including zero */
  npi = (flpt->p + 1) / 2;

  /* half-number of points, excluding zero */
  npx = flpt->p / 2;

  /* even transform */
  ar = flpt->ear;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_exp: insufficient working array\n");
    return;
  }

  /* copy even part of v to working array */

  /* treat zero specially */
  if (flpt->p % 2 == 1)
    w->v[ar->ncol + npi] = v->v[npi] * flpt->w->v[npi];

  /* copy other values */
  for (i=0; i< npx; i++)
    w->v[ar->ncol + i]
      = (v->v[i] + v->v[flpt->p - i - 1]) * flpt->w->v[i];

  /* do the even transform */
  fxt_sarld_exp(ar, w->v);

  /* copy results to u from working vector */
  for (i=0; i< (u->n + 1) / 2; i++)
    u->v[2 * i] = w->v[i];

  /* odd transform */
  if (flpt->n == 0)
    return;

  ar = flpt->oar;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_exp: insufficient working array\n");
    return;
  }

  /* copy odd part of v to working array */
  for (i=0; i< npx; i++)
    w->v[ar->ncol + i]
      = (v->v[i] - v->v[flpt->p - i - 1]) * flpt->w->v[i];

  /* do the odd transform */
  fxt_sarld_exp(ar, w->v);

  /* copy results to u from working array */
  for (i=0; i< u->n / 2; i++)
    u->v[2 * i + 1] = w->v[i];
}
