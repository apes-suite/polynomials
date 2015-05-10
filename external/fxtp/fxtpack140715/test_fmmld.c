#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "test_fxtpack.h"

#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"
#include "fxt_vecld.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing FMM of FXTPACK
*********************************************************************/

void test_fmmld(void) {
  int c;

  printf("- test_fmmld...\n");

  for (c = 0; c < 3; c ++) {
    fxt_vecld *x, *y, *ss, *st;
    fxt_matld *mat;
    fxt_fmmld *fmm;
    long n, m, i, j;
    double err;

    /* set different sizes */
    if (c == 0) {
      n = 40;  m = 60;
    } else if (c == 1) {
      n = 50;  m = 50;
    } else {
      n = 60;  m = 40;
    }

    /* allocate vectors */
    x = fxt_vecld_new(n);
    y = fxt_vecld_new(m);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* set random but sorted vector */
    x->v[0] = 0.0;
    for (i=1; i< n; i++)
      x->v[i] = x->v[i-1] + 2.0 * rand() / (n * (1.0 + RAND_MAX));

    y->v[0] = rand() / (m * (1.0 + RAND_MAX));
    for (j=1; j< m; j++)
      y->v[j] = y->v[j-1] + 2.0 * rand() / (m * (1.0 + RAND_MAX));

    /* create matrix */
    mat = fxt_matld_new(m, n);

    for (i=0; i< m; i++)
      for (j=0; j< n; j++)
	MENT(mat, i, j) = 1.0 / (y->v[i] - x->v[j]);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* make FMM */
    fmm = fxt_fmmld_new(mat, x, y, 1e-6);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    printf("- FMM %ldx%ld flop=%ld dir=%ld\n", fmm->ns, fmm->nt,
	   fxt_fmmld_evaluate_mflop(fmm), fmm->ns * fmm->nt);

    /* undo preprocessing... */
    fxt_fmmld_unpreproc(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* check estimation routines */
    test_fmmld_estimate(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* check preprocessing routines */
    test_fmmld_preproc(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* redo preprocessing */
    fxt_fmmld_preproc(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* check scaling */
    ss = fxt_vecld_new(fmm->ns);
    st = fxt_vecld_new(fmm->nt);

    /* random numbers [0.5, 1.5) */
    for (i=0; i< ss->n; i++)
      ss->v[i] = 0.5 + rand() / (1.0 + RAND_MAX);

    for (i=0; i< st->n; i++)
      st->v[i] = 0.5 + rand() / (1.0 + RAND_MAX);

    /* do scaling */
    fxt_fmmld_scalecol(fmm, ss);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_fmmld_scalerow(fmm, st);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    err = fxt_fmmld_error(fmm, POWERPREC);

    if (fmm->warn == 0 &&
	err > fmm->prec * 1.5 + 1e-16 * (fmm->ns + fmm->nt)) {
#if DISALLOWHIGHERROR
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld: error in scaling\n");
      return;
#else
      fxt_error_set(FXT_ERROR_WARN,
		    "test_fmmld: error in scaling\n");
#endif
    }

    fxt_vecld_del(st);
    fxt_vecld_del(ss);

    fxt_fmmld_del(fmm);

    fxt_matld_del(mat);
    fxt_vecld_del(y);
    fxt_vecld_del(x);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}
