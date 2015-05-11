#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test Givens LQ decomposition for matld
*********************************************************************/

void test_matld_glq(void) {
  int c;

  printf("- test_matld_glq...\n");

  for (c=0; c< 4; c++) {
    fxt_matld *a, *lq;
    fxt_vecl *piv;
    long n, m;
    long i, j;

    /* set sizes */
    if (c == 0) {
      n = 50;  m = 30;
    } else if (c == 1) {
      n = 50;  m = 20;
    } else if (c == 2) {
      n = 20;  m = 50;
    } else {
      n = 30;  m = 50;
    }

    a = fxt_matld_new(n, m);
    piv = fxt_vecl_new(n);

    for (i=0; i< n*m; i++)
      a->v[i] = (double) (rand() / (1.0 + RAND_MAX));

    lq = fxt_matld_clone(a, 0, n-1, 0, m-1);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    fxt_matld_glq(lq, piv);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    if (n < m) {
      fxt_matld *a0;

      fxt_matld_permrow(a, piv);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      fxt_matld_glq_mqinv(a, lq);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	for (j=n; j< m; j++)
	  if (fabs(MENT(a, i, j)) >= 1e-14) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_glq: error in M*Qinv\n");
	    return;
	  }

      a0 = fxt_matld_clone(a, 0, n-1, 0, n-1);
      fxt_matld_glq_mlinv(a0, lq);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	for (j=0; j< n; j++)
	  if (fabs(MENT(a0, i, j) - (i == j ? 1.0 : 0.0)) > 1e-14) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_glq: error in M*Linv\n");
	    return;
	  }
      
      fxt_matld_del(a0);
    } else {
      fxt_matld *a0, *a1, *l1;

      fxt_matld_permrow(a, piv);

      a0 = fxt_matld_clone(a, 0, m-1, 0, m-1);
      a1 = fxt_matld_clone(a, m, n-1, 0, m-1);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      fxt_matld_glq_mqinv(a0, lq);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      fxt_matld_glq_mlinv(a0, lq);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< m; i++)
	for (j=0; j< m; j++)
	  if (fabs(MENT(a0, i, j) - (i == j ? 1.0 : 0.0)) > 1e-14) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_glq: error in M*Ainv\n");
	    return;
	  }

      fxt_matld_glq_mqinv(a1, lq);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      l1 = fxt_matld_clone(lq, m, n-1, 0, m-1);
      fxt_matld_sub(a1, l1);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n-m; i++)
	for (j=0; j< m; j++)
	  if (fabs(MENT(a1, i, j)) > 1e-14) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_glq: error in M*Qinv\n");
	    return;
	  }

      fxt_matld_del(l1);
      fxt_matld_del(a1);
      fxt_matld_del(a0);

      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }

    fxt_matld_del(lq);
    fxt_vecl_del(piv);
    fxt_matld_del(a);

    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }
}
