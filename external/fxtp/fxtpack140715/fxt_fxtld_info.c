#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  evaluate resulting FXT
*********************************************************************/

#define MAX(x, y) ((x) > (y) ? (x) : (y))

static void part_evl(fxt_vecld *v, fxtld *fxt, fxt_vecld *u, long u0){
  fxt_vecld *w;			/* output vector */

  if (fxt->fmm != NULL) {
    w = fxt->ws;
    fxt_vecld_zero(w);
  } else
    w = v;

  if (fxt->de != NULL) {

    fxt_matld_addmatprtvec(w, fxt->de,
			   u, fxt->ui0 - u0, fxt->ui1 - u0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    part_evl(w, fxt->ch1, u, u0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    part_evl(w, fxt->ch0, u, u0);
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


static void part_exp(fxt_vecld *u, fxtld *fxt, fxt_vecld *v, long u0){
  fxt_vecld *w;

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

    fxt_matld_prtsettrmatvec(u, fxt->ui0 - u0, fxt->ui1 - u0,
			     fxt->de, w);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    part_exp(u, fxt->ch0, w, u0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    part_exp(u, fxt->ch1, w, u0);

  }
}


/*** evaluate error of FXT ***/
double fxtld_error(fxtld *fxt, fxt_matld *mat, double prec) {
  fxt_vecld *u, *uu, *v, *vv;
  double s, s0, ss;
  int c, cc;
  
  /* check null pointers */
  if (fxt == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_error: null pointer\n");
    return 0.0;
  }

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_error: irregal precision\n");
    return 0.0;
  }

  /* check size */
  if (mat->nrow == 0 || mat->ncol == 0)
    return 0.0;

  u  = fxt_vecld_new(mat->ncol);
  uu = fxt_vecld_new(mat->ncol);
  v  = fxt_vecld_new(mat->nrow);
  vv = fxt_vecld_new(mat->nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    fxt_vecld_rand(uu);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0.0;

    s = s0 = 0.0;

    for (c=0; ; c++) {
      fxt_vecld_normalize(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_setmatvec(v, mat, uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      if (fxt->par == NULL) {

	fxtld_evl(vv, fxt, uu);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0.0;

      } else {

	fxt_vecld_zero(vv);
	part_evl(vv, fxt, uu, fxt->ui0);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0.0;

      }

      fxt_vecld_sub(vv, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_settrmatvec(u, mat, vv);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      if (fxt->par == NULL) {

	fxtld_exp(uu, fxt, vv);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0.0;

      } else {

	part_exp(uu, fxt, vv, fxt->ui0);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0.0;

      }

      fxt_vecld_sub(uu, u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      s = fxt_vecld_norm2(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      if (s == 0.0 || (c > MINPOWERITER && (s - s0) / s <= prec))
	break;

      s0 = s;
    }

    s = MAX(s, s0);

    if (s == 0.0 || fabs(s - ss) / MAX(s, ss) <= prec)
      break;

    ss = MAX(s, ss);
  }

  s = MAX(s, ss);

  fxt_vecld_del(vv);
  fxt_vecld_del(v);
  fxt_vecld_del(uu);
  fxt_vecld_del(u);

  return sqrt(s);
}


/*** number of multiply operations ***/
long fxtld_evaluate_mflop(fxtld *fxt) {
  long flop;
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_evaluate_mflop: null pointer\n");
    return 0;
  }

  if (fxt->de != NULL) {

    flop = fxt_matld_dropflop(fxt->de, 0.0);

  } else {

    flop = fxtld_evaluate_mflop(fxt->ch0)
      + fxtld_evaluate_mflop(fxt->ch1);

  }

  if (fxt->fmm != NULL) {

    flop += fxt->ws->n + fxt_fmmld_evaluate_mflop(fxt->fmm);

  }

  return flop;
}


/*** number of temporal variables ***/
long fxtld_evaluate_work(fxtld *fxt) {
  long work = 0;
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_evaluate_work: null pointer\n");
    return 0;
  }

  if (fxt->fmm != NULL)
    work = fxt->ws->n + fxt->wt->n
      + fxt_fmmld_evaluate_work(fxt->fmm);

  if (fxt->ch0 != NULL)
    work += fxtld_evaluate_work(fxt->ch0)
      + fxtld_evaluate_work(fxt->ch1);

  return work;
}
