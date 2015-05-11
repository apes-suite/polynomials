#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "test_fxtpack.h"

#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"
#include "fxt_math.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"
#include "fxt_sarld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing FXT preprocessor
*********************************************************************/

/* check number of ins/outs */
static void check_n(fxt_lagld *lag) {
  long i;

  for (i=0; i< lag->nu_var; i++) {
    fxt_lagld_ent *e;  long n;
    fxt_lagld_var *var = &lag->valloc[i];

    n = var->io == 1 ? 1 : 0;
    for (e=var->ient; e != NULL; e=e->inext, n++);
    if (var->n_in != n) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fxtld: error in n_in\n");
      return;
    }

    n = var->io == 2 ? 1 : 0;
    for (e=var->oent; e != NULL; e=e->onext, n++);
    if (var->n_out != n) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fxtld: error in n_out\n");
      return;
    }

    if (var->io != 0 || var->n_in != 0 || var->n_out != 0)
      if (var->c != 0) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fxtld: uncleared counter\n");
	return;
      }
  }
}


static long check_fxt(fxt_matll *lg, fxt_vecll *xx, double eps) {
  fxtld *fxt;
  fxt_matld *dlg;
  double err;
  fxt_lagld *lag;
  long mflop0, mflop1, mflop2, mflop3;
  fxt_sarld *ar;
  
  fxt = fxtld_from_matll(lg, xx, eps);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  dlg = fxt_matll_to_matld(lg);
  mflop0 = lg->nrow * lg->ncol;

  err = fxtld_error(fxt, dlg, POWERPREC);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  if (err > eps) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: FXT error %e too large\n", err);
    return 0;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fxtld: FXT error %e too large\n", err);
#endif
  }

  mflop1 = fxtld_evaluate_mflop(fxt);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  if (mflop0 < mflop1) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: FXT operation increased\n");
    return 0;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fxtld: FXT operation increased\n");
#endif
  }

  lag = fxtld_get_lagld(fxt);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  check_n(lag);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  err = fxt_matld_lagld_error(lag, dlg, POWERPREC, 0.0);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  if (err > eps) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: LAG error %e too large\n", err);
    return 0;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fxtld: LAG error %e too large\n", err);
#endif
  }

  mflop2 = fxt_lagld_evaluate_mflop(lag);

  if (mflop1 < mflop2) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: LAG operation increased\n");
    return 0;
  }

  fxt_matld_lagld_opt(lag, dlg, eps);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  check_n(lag);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  err = fxt_matld_lagld_error(lag, dlg, POWERPREC, 0.0);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  if (err > eps) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: LAGOPT error %e too large\n", err);
    return 0;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fxtld: LAGOPT error %e too large\n", err);
#endif
  }

  mflop3 = fxt_lagld_evaluate_mflop(lag);
  if (mflop2 < mflop3) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: LAGOPT operation increased\n");
    return 0;
  }

  if (errmax < err) errmax = err;

  check_n(lag);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  ar = fxt_sarld_from_lagld(lag);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  err = fxt_sarld_error(ar, dlg, NULL, POWERPREC);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  if (err > eps) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fxtld: SAR error %e too large\n", err);
    return 0;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fxtld: SAR error %e too large\n", err);
#endif
  }

  fxt_sarld_del(ar);

  fxt_lagld_del(lag);
  fxt_matld_del(dlg);

  fxtld_del(fxt);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  return mflop3;
}

long test_fxtld(long p, long n, long m, double eps) {
  fxt_vecll *x, *w, *xx, *ww;
  fxt_matll *lg;
  long mflop = 0;

  printf("- test_fxtld(%ld %ld %ld %e)\n", p, n, m, eps);

  fxt_gauss_vecll(p, &x, &w);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;
  
  lg = fxt_legendre_matll(x, m, n, 0);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  xx = fxt_vecll_clone(x, 0, lg->nrow - 1);
  fxt_vecll_sq(xx);
  fxt_vecll_neg(xx);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
  fxt_vecll_scale(ww, 0, p/2 - 1, 2.0L);
  fxt_vecll_sqrt(ww);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  fxt_matll_scalerow(lg, ww);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  mflop += check_fxt(lg, xx, eps);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  fxt_vecll_del(ww);
  fxt_vecll_del(xx);
  fxt_matll_del(lg);

  if (m != n) {
    lg = fxt_legendre_matll(x, m, n, 1);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return 0;

    xx = fxt_vecll_clone(x, 0, lg->nrow - 1);
    fxt_vecll_sq(xx);
    fxt_vecll_neg(xx);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return 0;

    ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
    fxt_vecll_scale(ww, 0, p/2 - 1, 2.0L);
    fxt_vecll_sqrt(ww);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return 0;

    fxt_matll_scalerow(lg, ww);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return 0;

    mflop += check_fxt(lg, xx, eps);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return 0;

    fxt_vecll_del(ww);
    fxt_vecll_del(xx);
    fxt_matll_del(lg);
  }

  fxt_vecll_del(w);
  fxt_vecll_del(x);

  if (fxt_error_level() > FXT_ERROR_WARN)
    return 0;

  return mflop;
}
