#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "test_fxtpack.h"

#include "fxt_flptld.h"
#include "fxt_flptld_loc.h"
#include "fxt_math.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test loaded Fast Legendre Polynomial Transform
*********************************************************************/

static void check_legendre(int odd, fxt_vecll *x, fxt_vecll *w,
			   fxt_flptld *flpt) {
  long n;
  fxt_matll *lg;
  fxt_matld *dlg;
  fxt_vecll *ww;
  fxt_vecld *dw;
  fxt_sarld *ar;
  double err;

  n = flpt->n;
  ar = (odd ? flpt->oar : flpt->ear);

  if (odd && n == 0)
    return;

  lg = fxt_legendre_matll(x, 0, n, odd);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecll_scale(ww, 0, w->n / 2 - 1, 2.0L);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  dw = fxt_vecll_to_vecld(ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  dlg = fxt_matll_to_matld(lg);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  err = fxt_sarld_error(ar, dlg, dw, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (err > flpt->prec) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_flptld: error %e too large\n", err);
    return;
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "test_flptld: error %e too large\n", err);
#endif
  } 

  fxt_matld_del(dlg);
  fxt_vecld_del(dw);
  fxt_vecll_del(ww);
  fxt_matll_del(lg);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}


void test_flptld(char *fname) {
  fxt_flptld *flpt;
  fxt_vecll *x, *w;
  fxt_vecld *xx, *ww;

  printf("- testing %s\n", fname);

  flpt = fxt_flptld_load(fname);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_gauss_vecll(flpt->p, &x, &w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  xx = fxt_vecll_to_vecld(x);
  ww = fxt_vecll_to_vecld(w);

  fxt_vecld_sub(xx, flpt->x);
  fxt_vecld_sub(ww, flpt->w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (fxt_vecld_norm2(xx) / fxt_vecld_norm2(flpt->x) > 1e-15) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_flptld: Gauss point mismatch\n");
    return;
  }

  if (fxt_vecld_norm2(ww) / fxt_vecld_norm2(flpt->w) > 1e-15) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_flptld: Gauss weight mismatch\n");
    return;
  }

  fxt_vecld_del(ww);
  fxt_vecld_del(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  check_legendre(0, x, w, flpt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  check_legendre(1, x, w, flpt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecll_del(w);
  fxt_vecll_del(x);

  fxt_flptld_del(flpt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}
