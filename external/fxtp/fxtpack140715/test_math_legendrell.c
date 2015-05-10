#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_math.h"
#include "fxt_matll.h"
#include "fxt_vecll.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test associated Legendre function generation
*********************************************************************/

void test_math_legendrell(void) {
  long n, m, odd;

  printf("- test_math_legendrell...\n");

  for (n = 50; n < 60; n++) {
    fxt_vecll *x, *w;

    fxt_gauss_vecll(n, &x, &w);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (m = 0; m < n; m ++) {
      for (odd = 0; odd <= 1; odd ++) {
	fxt_matll *lg;
	fxt_vecll *ww;
	fxt_matld *ld;
	fxt_vecld *c0, *c1;
	long i, j;
	double dd;

	if (m == n-1 && odd)
	  continue;

	lg = fxt_legendre_matll(x, m, n-1, odd);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	ww = fxt_vecll_clone(w, 0, lg->nrow - 1);
	fxt_vecll_scale(ww, 0, n/2 - 1, 2.0L);
	fxt_vecll_sqrt(ww);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	fxt_matll_scalerow(lg, ww);
	fxt_vecll_del(ww);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	ld = fxt_matll_to_matld(lg);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	c0 = fxt_vecld_new(ld->nrow);
	c1 = fxt_vecld_new(ld->nrow);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< ld->ncol; i++) {
	  fxt_matld_getcol(c0, ld, i);
	  if (fxt_error_level() > FXT_ERROR_WARN)
	    return;

	  for (j=0; j< ld->ncol; j++) {
	    fxt_matld_getcol(c1, ld, j);
	    if (fxt_error_level() > FXT_ERROR_WARN)
	      return;

	    dd = (i == j ? 1.0 : 0.0) - fxt_vecld_dot(c0, c1);

	    if (fabs(dd) > 2e-14 * ld->ncol) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_math_legendrell:"
			    " orthogonality error\n");
	      return;
	    }
	  }
	}

	fxt_vecld_del(c1);
	fxt_vecld_del(c0);
	fxt_matld_del(ld);

	fxt_matll_del(lg);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
      }
    }

    fxt_vecll_del(w);
    fxt_vecll_del(x);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }
}
