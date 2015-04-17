#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  evaluate resulting FMM
*********************************************************************/

double fxt_fmmld_error(fxt_fmmld *fmm, double prec) {
  fxt_vecld *u, *uu, *v, *vv;
  double s, s0, ss;
  int c, cc;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_error: null pointer\n");
    return 0.0;
  }

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_error: irregal precision\n");
    return 0.0;
  }

  /* check size */
  if (fmm->ns == 0 || fmm->nt == 0)
    return 0.0;

  /* allocate vectors */
  u  = fxt_vecld_new(fmm->ns);
  uu = fxt_vecld_new(fmm->ns);
  v  = fxt_vecld_new(fmm->nt);
  vv = fxt_vecld_new(fmm->nt);

  if (u == NULL || uu == NULL || v == NULL || vv == NULL) {
    fxt_error_raise();
    return 0.0;
  }

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    do {

      /* initialize by random vector */
      fxt_vecld_rand(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

    } while (fxt_vecld_norm2(uu) == 0.0);

    /* initialize norm */
    s = s0 = 0.0;

    /* power method iteration */
    for (c = 0; ; c ++) {

      /* normalize input */
      fxt_vecld_normalize(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* v = A uu */
      fxt_matld_setmatvec(v, fmm->mat, uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* vv = B uu */
      fxt_fmmld_evl(vv, fmm, uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* vv = vv - v = (B - A) uu */
      fxt_vecld_sub(vv, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* u = A^T vv */
      fxt_matld_settrmatvec(u, fmm->mat, vv);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* uu = B^T vv */
      fxt_fmmld_exp(uu, fmm, vv);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* uu = uu - u = (B - A)^T (B - A) uu */
      fxt_vecld_sub(uu, u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* get the norm */
      s = fxt_vecld_norm2(uu);
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

  /* deallocate temporals */
  fxt_vecld_del(vv);
  fxt_vecld_del(v);
  fxt_vecld_del(uu);
  fxt_vecld_del(u);

  /* return norm */
  return sqrt(s);
}


static double ksum;
static long nk;
static int kmax;

static void eval_reg(fmmld_reg *reg) {
  if (reg->type & FMM_RTMUL) {
    if (reg->km <= (reg->ns + 1) / 2) {
      nk ++;
      ksum += reg->km * reg->km;
    }

    if (kmax < reg->km)
      kmax = reg->km;
  }

  if (reg->type & FMM_RTLOC) {
    if (reg->kl <= (reg->nt + 1) / 2) {
      nk ++;
      ksum += reg->kl * reg->kl;
    }

    if (kmax < reg->kl)
      kmax = reg->kl;
  }

  if (reg->ch0 != NULL) {
    eval_reg(reg->ch0);
    eval_reg(reg->ch1);
  }
}

void fxt_fmmld_evaluate_kmax(fxt_fmmld *fmm) {
  int i;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_evaluate_kmax: null pointer\n");
    return;
  }

  ksum = 0.0;
  nk = 0;
  kmax = 0;

  for (i=0; i< fmm->n_roots; i++)
    eval_reg(&fmm->roots[i]);

  if (nk == 0)
    fmm->kest = -1;
  else
    fmm->kest = (int) (sqrt(ksum / nk) + 0.5);

  fmm->kmax = kmax;

  fmm->min_np = fmm->kest * 2;
  if (fmm->min_np < 2)
    fmm->min_np = 2;
}


static long evaluate_mflop(fmmld_reg *reg) {
  fmmld_cell *cell;
  long flop = 0;

  if (reg->pm != NULL)
    flop += fxt_matld_dropflop(reg->pm, 0.0);

  if (reg->lp != NULL)
    flop += fxt_matld_dropflop(reg->lp, 0.0);

  if (reg->ll != NULL)
    flop += fxt_matld_dropflop(reg->ll, 0.0);

  if (reg->mm != NULL)
    flop += fxt_matld_dropflop(reg->mm, 0.0);

  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next)
    flop += fxt_matld_dropflop(cell->m, 0.0);

  if (reg->ch0 != NULL) {
    flop += evaluate_mflop(reg->ch0);
    flop += evaluate_mflop(reg->ch1);
  }

  return flop;
}


long fxt_fmmld_evaluate_mflop(fxt_fmmld *fmm) {
  int i;
  long flop = 0;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_evaluate_mflop: null pointer\n");
    return 0;
  }

  /* sum-up mflops */
  for (i=0; i< fmm->n_roots; i++)
    flop += evaluate_mflop(&fmm->roots[i]);

  return flop;
}


static long evaluate_work(fmmld_reg *reg) {
  long work = 0;

  if (reg->type & FMM_RTMUL)
    work += reg->km;

  if (reg->type & FMM_RTLOC)
    work += reg->kl;

  if (reg->ch0 != NULL)
    work += evaluate_work(reg->ch0) + evaluate_work(reg->ch1);
  
  return work;
}


long fxt_fmmld_evaluate_work(fxt_fmmld *fmm) {
  int i;
  long work = 0;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_evaluate_work: null pointer\n");
    return 0;
  }

  for (i=0; i< fmm->n_roots; i++)
    work += evaluate_work(&fmm->roots[i]);

  return work;
}
