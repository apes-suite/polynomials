#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_flptld.h"
#include "fxt_flptld_loc.h"

#include "fxt_math.h"
#include "fxt_fxtld.h"
#include "fxt_lagld.h"
#include "fxt_sarld.h"

#include "fxt_matll.h"
#include "fxt_vecll.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generate preprocessed file for Fast Legendre Polynomial Transform
*********************************************************************/

static double res_prec;
static long res_mflop;

fxt_sarld* create_leg_half(int odd, fxt_vecll *x, fxt_vecll *w,
			   long n, double prec) {
  fxt_matll *lg;		/* the Legendre matrix */
  fxt_vecll *xx, *ww;		/* halved vectors */
  fxt_matld *dlg;		/* double-valued matrix */
  fxtld *fxt;			/* FXT structure */
  fxt_lagld *lag;		/* linear algorithm graph */
  fxt_sarld *ar;		/* simple array representation */
  fxt_vecld *dw;
  double err;

  /* skip odd part for the 'last' order */
  if (odd && 0 == n)
    return NULL;

  /* get Legendre function matrix */
  lg = fxt_legendre_matll(x, 0, n, odd);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* get half of the Gaussian points */
  xx = fxt_vecll_clone(x, 0, lg->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* functions are polynomials of x**2 */
  fxt_vecll_sq(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* negate for upward sorting */
  fxt_vecll_neg(xx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* get the halved Gaussian weight vector */
  ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* double the value of duplicated points */
  fxt_vecll_scale(ww, 0, w->n / 2 - 1, 2.0L);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* square-roots, being multiplied twice */
  fxt_vecll_sqrt(ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* scale the Legendre matrix */
  fxt_matll_scalerow(lg, ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* create fxt: the main part */
  fxt = fxtld_from_matll(lg, xx, prec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* get the double-valued matrix */
  dlg = fxt_matll_to_matld(lg);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* convert to linear algorithm graph */
  lag = fxtld_get_lagld(fxt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* optimize as for the precision */
  fxt_matld_lagld_opt(lag, dlg, prec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* check error */
  err = fxt_matld_lagld_error(lag, dlg, POWERPREC, 0.0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (err > prec) {
#if DISALLOWHIGHERROR
    if (err > prec + 1e-16 * lg->nrow) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_flptld_preproc: precision unattained\n");
      return NULL;
    } else
#endif
      fxt_error_set(FXT_ERROR_WARN,
		    "fxt_flptld_preproc: precision unattained\n");
  }

  res_prec = err;
  res_mflop = fxt_lagld_evaluate_mflop(lag);

  /* convert to simple array representation */
  ar = fxt_sarld_from_lagld(lag);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* descale sarld */
  dw = fxt_vecll_to_vecld(ww);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_einv(dw);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_sarld_scalerow(ar, dw);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_del(dw);


  /* deallocate everything */
  fxt_lagld_del(lag);
  fxt_matld_del(dlg);
  fxtld_del(fxt);
  fxt_vecll_del(ww);
  fxt_vecll_del(xx);
  fxt_matll_del(lg);

  return ar;
}

static void create_legendre(int odd, fxt_vecll *x, fxt_vecll *w,
			    long n, double prec, FILE *fout) {
  fxt_sarld *ar;		/* simple array representation */
  double err;


  ar = create_leg_half(odd, x, w, n, prec);

  /* save the results */
  fxt_sarld_save(fout, ar);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}

// Initialize a flptld data structure
fxt_flptld* fxt_flptld_init(long p, long n, double prec) {
  fxt_flptld *flpt;
  fxt_vecll *x, *w;

  /* check input */
  if (p < n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_init: irregal degree n\n");
    return NULL;
  }

  /* allocate body */
  flpt = (fxt_flptld*) malloc(sizeof(fxt_flptld));
  if (flpt == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_flptld_load: allocation failed\n");
    return NULL;
  }

  flpt->p = p;
  flpt->n = n;
  flpt->prec = prec;

  fxt_gauss_vecll(p, &x, &w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* Gaussian points */
  flpt->x = fxt_vecll_to_vecld(x);

  /* Gaussian weights */
  flpt->w = fxt_vecll_to_vecld(w);

  /* create even half */
  flpt->ear = create_leg_half(0, x, w, n, prec);

  /* create odd half */
  flpt->oar = create_leg_half(1, x, w, n, prec);

  /* deallocate vectors */
  fxt_vecll_del(w);
  fxt_vecll_del(x);

  return flpt;
}

void fxt_flptld_preproc(long p, long n, double prec, char *fname) {
  fxt_vecll *x, *w;
  fxt_vecld *v;
  long mflop0 = 0, dflop0 = 0;
  double errmax0 = 0.0;
  long size[2];
  FILE *fout;

  /* check null pointers */
  if (fname == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_preproc: null pointer\n");
    return;
  }

  /* check input */
  if (p < n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fxtld_preproc: irregal degree n\n");
    return;
  }

  /* open file */
  fout = fopen(fname, "wb");
  if (fout == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fxtld_preproc: cannot open %s\n", fname);
    return;
  }

  /* write header */
  fxt_file_writestr(fout, fxt_flptld_header);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* write size */
  size[0] = p;
  size[1] = n;
  fxt_file_writelongs(fout, size, 2);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* write precision */
  fxt_file_writedoubles(fout, &prec, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* create Gaussian points */
  fxt_gauss_vecll(p, &x, &w);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* save Gaussian points */
  v = fxt_vecll_to_vecld(x);

  fxt_vecld_save(fout, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_del(v);

  /* save Gaussian weights */
  v = fxt_vecll_to_vecld(w);

  fxt_vecld_save(fout, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_vecld_del(v);

  printf("- Legendre(%ld) %ld points\n", n, p);
  fflush(stdout);

  /* create even half */
  create_legendre(0, x, w, n, prec, fout);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  mflop0 += res_mflop;
  if (errmax0 < res_prec)
    errmax0 = res_prec;

  /* create odd half */
  create_legendre(1, x, w, n, prec, fout);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  mflop0 += res_mflop;
  if (errmax0 < res_prec)
    errmax0 = res_prec;

  /* total results */
  dflop0 = (n + 1) * p / 2;

  printf("- error=%9.2e (req_prec=%9.2e) speedup=%4.2f\n",
	 errmax0, prec, (double) dflop0 / (double) mflop0);

  /* finalize file */
  fxt_file_writestr(fout, fxt_flptld_header);

  /* deallocate vectors */
  fxt_vecll_del(w);
  fxt_vecll_del(x);

  /* close file */
  fclose(fout);
}
