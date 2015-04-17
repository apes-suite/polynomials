#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_math.h"
#include "fxt_vecll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test Gaussian nodes
*********************************************************************/

void test_math_gaussll(void) {
  fxt_vecll *xv, *wv;
  long n, i;
  long double v, eps;

  printf("- test_math_gaussll...\n");

  /* compute epsilon for long double */
  for (eps = 1.0L; 1.0L + eps != 1.0L; eps *= 0.5);

  for (n = 10; n < 100; n += 13) {
    fxt_gauss_vecll(n, &xv, &wv);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (i=0; i< n; i++) {
      if (-1.0L >= xv->v[i] || xv->v[i] >= 1.0L) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_math_gaussll: node range error\n");
	return;
      }

      if (0.0L >= wv->v[i] || wv->v[i] > 1.0L) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_math_gaussll: weight range error\n");
	return;
      }

      if (i > 0 && xv->v[i-1] >= xv->v[i]) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_math_gaussll: node order error\n");
	return;
      }
    }

    v = 0.0L;
    for (i=0; i< n; i++)
      v += wv->v[i];

    if (fabs((double)(v - 1.0L)) >= 30 * n * eps) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_math_gaussll: 0th order error\n");
      return;
    }

    v = 0.0L;
    for (i=0; i< n; i++)
      v += xv->v[i] * wv->v[i];

    if (fabs((double)v) >= 30 * n * eps) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_math_gaussll: 1st order error\n");
      return;
    }

    fxt_vecll_del(wv);
    fxt_vecll_del(xv);
  }
}
