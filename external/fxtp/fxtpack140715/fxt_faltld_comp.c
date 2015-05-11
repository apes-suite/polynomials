#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_faltld.h"
#include "fxt_faltld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  Fast Associated Legendre Transform
*********************************************************************/

void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
		    fxt_vecld *u, fxt_vecld *w) {
  long npi, npx, i;
  fxt_sarld *ar;
  double t;

  /* check pointers */
  if (v == NULL || falt == NULL || u == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_evl: null pointer\n");
    return;
  }

  /* check size */
  if (falt->p != v->n || falt->n - m + 1 < u->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_evl: size mismatch\n");
    return;
  }

  /* check order m */
  if (m < 0 || falt->n < m || falt->alt[m] == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_evl: irregal order m\n");
    return;
  }

  /* half-number of points, including equator */
  npi = (falt->p + 1) / 2;

  /* half-number of points, excluding equator */
  npx = falt->p / 2;

  /*** even transform ***/
  ar = falt->alt[m]->ear;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_evl: insufficient working array\n");
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

  /* treat equator specially */
  if (falt->p % 2 == 1)
    v->v[npi] = w->v[ar->ncol + npi];

  /* copy other values */
  for (i=0; i< npx; i++) {
    t = w->v[ar->ncol + i];
    v->v[i] = t;
    v->v[falt->p - i - 1] = t;
  }

  /*=== odd transform ===*/
  if (m == falt->n)
    return;

  ar = falt->alt[m]->oar;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_evl: insufficient working array\n");
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
    v->v[falt->p - i - 1] -= t;
  }
}

void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
		    fxt_vecld *v, fxt_vecld *w) {
  long npi, npx, i;
  fxt_sarld *ar;

  /* check pointers */
  if (v == NULL || falt == NULL || u == NULL || w == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_exp: null pointer\n");
    return;
  }

  /* check size */
  if (falt->p != v->n || falt->n - m + 1 < u->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_exp: size mismatch\n");
    return;
  }

  /* check order m */
  if (m < 0 || falt->n < m || falt->alt[m] == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_exp: irregal order m\n");
    return;
  }

  /* half-number of points, including equator */
  npi = (falt->p + 1) / 2;

  /* half-number of points, excluding equator */
  npx = falt->p / 2;

  /* even transform */
  ar = falt->alt[m]->ear;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_exp: insufficient working array\n");
    return;
  }

  /* copy even part of v to working array */

  /* treat equator specially */
  if (falt->p % 2 == 1)
    w->v[ar->ncol + npi] = v->v[npi] * falt->w->v[npi];

  /* copy other values */
  for (i=0; i< npx; i++)
    w->v[ar->ncol + i]
      = (v->v[i] + v->v[falt->p - i - 1]) * falt->w->v[i];

  /* do the even transform */
  fxt_sarld_exp(ar, w->v);

  /* copy results to u from working vector */
  for (i=0; i< (u->n + 1) / 2; i++)
    u->v[2 * i] = w->v[i];

  /* odd transform */
  if (m == falt->n)
    return;

  ar = falt->alt[m]->oar;

  /* check working vector */
  if (w->n < ar->ncol + ar->nrow + ar->ntmp) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_exp: insufficient working array\n");
    return;
  }

  /* copy odd part of v to working array */
  for (i=0; i< npx; i++)
    w->v[ar->ncol + i]
      = (v->v[i] - v->v[falt->p - i - 1]) * falt->w->v[i];

  /* do the odd transform */
  fxt_sarld_exp(ar, w->v);

  /* copy results to u from working array */
  for (i=0; i< u->n / 2; i++)
    u->v[2 * i + 1] = w->v[i];
}
