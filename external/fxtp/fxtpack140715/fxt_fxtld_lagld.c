#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"
#include "fxt_fmmld.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  FXT to LAG data structure
*********************************************************************/


static long eval_mflop(fxtld *fxt) {
  long mflop;

  if (fxt->de != NULL)
    mflop = fxt_matld_dropflop(fxt->de, 0.0);
  else
    mflop = eval_mflop(fxt->ch0) + eval_mflop(fxt->ch1);

  if (fxt->fmm != NULL)
    mflop += fxt->ws->n + fxt->wt->n
      + fxt_fmmld_evaluate_mflop(fxt->fmm);

  return mflop;
}


/*** create linear algorithm graph ***/
fxt_lagld *fxtld_get_lagld(fxtld *fxt) {
  fxt_lagld *lag;
  long mflop, work;
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_get_lagld: null pointer\n");
    return NULL;
  }
  
  mflop = eval_mflop(fxt);
  work = fxtld_evaluate_work(fxt);

  /* allocate data structure */
  lag = fxt_lagld_new(fxt->mat->ncol, fxt->mat->nrow,
		      work, mflop);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* setup linear algrithm graph */
  fxtld_lagld(lag, fxt_lagld_ovec(lag), fxt, fxt_lagld_ivec(lag));
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* normalize linear algorithm graph */
  fxt_lagld_normalize(lag);

  return lag;
}


/*** set linear algorithm graph ***/
void fxtld_lagld(fxt_lagld *lag, fxt_lagld_vec *v,
		     fxtld *fxt, fxt_lagld_vec *u) {
  fxt_lagld_vec *w;

  if (fxt->fmm != NULL) {

    w = fxt_lagld_newvec(lag, fxt->ws->n);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else
    w = v;

  if (fxt->de != NULL) {

    fxt_matld_lagld(lag, w, 0, fxt->de, u, fxt->ui0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    fxtld_lagld(lag, w, fxt->ch0, u);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxtld_lagld(lag, w, fxt->ch1, u);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  }

  if (fxt->fmm != NULL) {
    fxt_lagld_vec *wt;
    long i, ns, nt;

    ns = fxt->ws->n;
    nt = fxt->wt->n;

    wt = fxt_lagld_newvec(lag, nt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    for (i=0; i< ns; i++) {
      fxt_lagld_add(lag, &v->v[fxt->idx->v[i]], 1.0, &w->v[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    fxt_fmmld_lagld(lag, wt, fxt->fmm, w);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    for (i=0; i< nt; i++) {
      fxt_lagld_add(lag, &v->v[fxt->idx->v[i + ns]], 1.0, &wt->v[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    free(wt);
    free(w);
  }
}
