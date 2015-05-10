#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fxtld_loc.h"
#include "fxt_fxtld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  evaluation and expansion of FXT
*********************************************************************/

/*** evaluate FXT ***/
void fxtld_evl(fxt_vecld *v, fxtld *fxt, fxt_vecld *u) {
  fxt_vecld *w;
  
  /* check null pointers */
  if (v == NULL || fxt == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_evl: null pointer\n");
    return;
  }

  if (fxt->par == NULL)
    fxt_vecld_zero(v);

  if (fxt->fmm != NULL) {
    w = fxt->ws;
    fxt_vecld_zero(w);
  } else
    w = v;

  if (fxt->de != NULL) {

    fxt_matld_addmatprtvec(w, fxt->de, u, fxt->ui0, fxt->ui1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    fxtld_evl(w, fxt->ch0, u);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxtld_evl(w, fxt->ch1, u);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  }

  if (fxt->fmm != NULL) {
    long i, ns, nt;

    ns = fxt->ws->n;
    nt = fxt->wt->n;

    for (i=0; i< ns; i++)
      v->v[fxt->idx->v[i]] += fxt->ws->v[i];

    fxt_fmmld_evl(fxt->wt, fxt->fmm, fxt->ws);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    for (i=0; i< nt; i++)
      v->v[fxt->idx->v[i + ns]] += fxt->wt->v[i];
  }
}

/*** expand (transposed evaluate) FXT ***/
void fxtld_exp(fxt_vecld *u, fxtld *fxt, fxt_vecld *v) {
  fxt_vecld *w;
  
  /* check null pointers */
  if (u == NULL || fxt == NULL || v == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_exp: null pointer\n");
    return;
  }

  if (fxt->fmm != NULL) {
    long i, ns, nt;

    ns = fxt->ws->n;
    nt = fxt->wt->n;

    for (i=0; i< nt; i++)
      fxt->wt->v[i] = v->v[fxt->idx->v[i + ns]];

    fxt_fmmld_exp(fxt->ws, fxt->fmm, fxt->wt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    for (i=0; i< ns; i++)
      fxt->ws->v[i] += v->v[fxt->idx->v[i]];

    w = fxt->ws;
  } else
    w = v;

  if (fxt->de != NULL) {

    fxt_matld_prtsettrmatvec(u, fxt->ui0, fxt->ui1, fxt->de, w);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    fxtld_exp(u, fxt->ch0, w);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxtld_exp(u, fxt->ch1, w);
  }
}
