#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  A-matrix stuff for FXT
*********************************************************************/

/*** the size of output vector ***/
static long get_m(fxtld *fxt) {
  long m;

  /* no info for NULL data */
  if (fxt == NULL)
    return 0;

  /* get data from parent */
  m = get_m(fxt->par);

  /* size of the first FMM */
  if (m == 0 && fxt->fmm != NULL)
    return fxt->idx->n;

  /* otherwise, data from parent */
  return m;
}


/*** get the output vector ***/
static fxt_vecld* get_v(fxtld *fxt, fxt_vecld *v) {

  /* go over the root */
  if (fxt == NULL)
    return v;

  /* interpolation temporal */
  if (fxt->fmm != NULL)
    return fxt->ws;

  /* otherwise, resort to parent */
  return get_v(fxt->par, v);
}


/*** evaluation by reverse recursion ***/
static void rev_evl(fxtld *fxt, fxt_vecld *v0) {

  /* get the output vector for this level */
  fxt_vecld *v = get_v(fxt->par, v0);

  if (fxt->fmm != NULL) {
    /* interpolation */

    long i, ns, nt, n;

    ns = fxt->ws->n;
    nt = fxt->wt->n;
    n = fxt->idx->n;

    /* ws must be already set */
    for (i=0; i< ns; i++)
      v->v[fxt->idx->v[i]] = fxt->ws->v[i];

    /* interpolate to wt */
    fxt_fmmld_evl(fxt->wt, fxt->fmm, fxt->ws);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* copy interpolated result */
    for (i=0; i< nt; i++)
      v->v[fxt->idx->v[i + ns]] = fxt->wt->v[i];

    /* clear dropped elements */
    for (i= ns + nt; i < n; i++)
      v->v[fxt->idx->v[i]] = 0.0;
  }

  /* reverse recursion to parent */
  if (fxt->par != NULL)
    rev_evl(fxt->par, v0);
}


/*** expansion by reverse recursion ***/
static fxt_vecld* rev_exp(fxtld *fxt, fxt_vecld *v) {

  /* reverse recursion to parent */
  if (fxt->par != NULL) {
    v = rev_exp(fxt->par, v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;
  }

  if (fxt->fmm != NULL) {
    long i, ns, nt;

    ns = fxt->ws->n;
    nt = fxt->wt->n;

    /* copy it to wt */
    for (i=0; i< nt; i++)
      fxt->wt->v[i] = v->v[fxt->idx->v[i + ns]];

    /* interpolation */
    fxt_fmmld_exp(fxt->ws, fxt->fmm, fxt->wt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* store results */
    for (i=0; i< ns; i++)
      fxt->ws->v[i] += v->v[fxt->idx->v[i]];

    /* the result of this level */
    v = fxt->ws;
  }

  /* return result */
  return v;
}


/*** get A-scale ***/
fxt_vecld *fxtld_get_ascale(fxtld *fxt) {
  long m;  fxt_vecld *a, *u, *v;  long i;
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_get_ascale: null pointer\n");
    return NULL;
  }

  /* get the output size */
  m = get_m(fxt);

  /* no output: this must be an error */
  if (m == 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_get_ascale: no interpolation\n");
    return NULL;
  }

  /* try output vector */
  u = get_v(fxt, NULL);

  /* allocate result vector */
  a = fxt_vecld_new((u == NULL) ? m : u->n);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* allocate output vector */
  v = fxt_vecld_new(m);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* get next output vector */
  u = get_v(fxt, v);

  /* simply call with unit vector */
  for (i=0; i< u->n; i++) {

    /* set unit vector */
    fxt_vecld_unit(u, i);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* clear result */
    fxt_vecld_zero(v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* evaluate */
    rev_evl(fxt, v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* get scale value */
    a->v[i] = fxt_vecld_norm2(v);
  }

  /* deallocate temporal */
  fxt_vecld_del(v);

  /* return A-scale vector */
  return a;
}


/*** get A-norm ***/
double fxtld_get_anorm(fxtld *fxt, double prec) {
  long m;  int c, cc;
  fxt_vecld *a, *u, *v, *w;
  double s, s0, ss;
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_get_anorm: null pointer\n");
    return 0.0;
  }

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_get_anorm: irregal precision\n");
    return 0.0;
  }

  /* get the output size */
  m = get_m(fxt);

  /* check for no computation (perhaps fxt is root) */
  if (m == 0)
    return 1.0;

  /* allocate vector */
  v = fxt_vecld_new(m);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  /* get the A-scale vector */
  a = fxtld_get_ascale(fxt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  /* invert for descaling */
  fxt_vecld_einv(a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  /* get the first level result vector */
  u = get_v(fxt, v);

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    do {

      /* initialize input as random */
      fxt_vecld_rand(u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

    } while (fxt_vecld_norm2(u) == 0.0);

    /* initialize norms */
    s = s0 = 0.0;

    /* power method iteration */
    for (c=0; ; c++) {

      /* normalize input */
      fxt_vecld_normalize(u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* scale u */
      fxt_vecld_emul(u, a);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* clear resulting vector */
      fxt_vecld_zero(v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* evaluate */
      rev_evl(fxt, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* expansion */
      w = rev_exp(fxt, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* scale w */
      fxt_vecld_emul(w, a);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* get the norm */
      s = fxt_vecld_norm2(w);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* check convergence */
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

  /* deallocate vectors */
  fxt_vecld_del(a);
  fxt_vecld_del(v);

  /* return result */
  return sqrt(s);
}
