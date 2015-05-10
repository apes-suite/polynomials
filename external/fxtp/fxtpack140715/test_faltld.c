#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "test_fxtpack.h"

#include "fxt_faltld.h"
#include "fxt_faltld_loc.h"
#include "fxt_math.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing Fast Associated Legendre Transform
*********************************************************************/

static double check_legendre(int odd, long n, long m,
			     fxt_vecll *x, fxt_vecll *w,
			     fxt_sarld *ar) {
  double err;
  fxt_matll *lg;
  fxt_vecll *ww;
  fxt_vecld *dw;
  fxt_matld *dlg;

  if (odd && n == m)
    return 0.0;

  lg = fxt_legendre_matll(x, m, n, odd);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
  fxt_vecll_scale(ww, 0, w->n / 2 - 1, 2.0L);
  dw = fxt_vecll_to_vecld(ww);

  dlg = fxt_matll_to_matld(lg);

  err = fxt_sarld_error(ar, dlg, dw, POWERPREC);

  fxt_matld_del(dlg);
  fxt_vecld_del(dw);
  fxt_vecll_del(ww);
  fxt_matll_del(lg);

  return err;
}

void test_faltld(char *fname) {
  fxt_faltld *falt;
  fxt_vecll *x, *w;
  fxt_vecld *xx, *ww;
  long m;
  double emax;

  printf("- testing %s\n", fname);

  falt = fxt_faltld_load(fname);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_gauss_vecll(falt->p, &x, &w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  xx = fxt_vecll_to_vecld(x);
  ww = fxt_vecll_to_vecld(w);

  fxt_vecld_sub(xx, falt->x);
  fxt_vecld_sub(ww, falt->w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (fxt_vecld_norm2(xx) / fxt_vecld_norm2(falt->x) > 1e-15) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_faltld: Gauss point mismatch\n");
    return;
  }

  if (fxt_vecld_norm2(ww) / fxt_vecld_norm2(falt->w) > 1e-15) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_faltld: Gauss weight mismatch\n");
    return;
  }

  fxt_vecld_del(ww);
  fxt_vecld_del(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  emax = 0.0;
  
  for (m=0; m<= falt->n; m++)
    if (falt->alt[m] != NULL) {
      double err0, err1;

      err0 = check_legendre(0, falt->n, m, x, w,
			    falt->alt[m]->ear);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      err1 = check_legendre(1, falt->n, m, x, w,
			    falt->alt[m]->oar);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      if (err0 < err1)
	err0 = err1;

      if (err0 > falt->prec) {
#if DISALLOWHIGHERROR
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_faltld: m=%ld err=%9.2e\n", m, err0);
	return;
#else
	fxt_error_set(FXT_ERROR_WARN,
		      "test_faltld: m=%ld err=%9.2e\n", m, err0);
#endif
      }

      if (emax < err0)
	emax = err0;
    }

  fxt_vecll_del(w);
  fxt_vecll_del(x);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  test_faltld_comp(falt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_faltld_del(falt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}
