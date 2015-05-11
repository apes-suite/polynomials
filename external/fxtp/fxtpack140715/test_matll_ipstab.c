#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matll.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test stable interpolation for matll
*********************************************************************/

void test_matll_ipstab(void) {
  long n, m, i, j;
  fxt_matll *a, *a0, *a1, *r, *cc;
  fxt_matld *c;
  fxt_vecl *p;
  int ccc;

  printf("- test_matll_ipstab...\n");

  for (ccc = 0; ccc < 2; ccc ++) {
    /* set size */
    if (ccc == 0) {
      n = 50;  m = 20;
    } else {
      n = 50;  m = 30;
    }

    a = fxt_matll_new(n, m);
    for (i=0; i< n*m; i++)
      a->v[i] = (long double) (rand() / (1.0 + RAND_MAX));

    p = fxt_vecl_new(n);

    c = fxt_matld_new(n - m, m);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    fxt_matll_ipstab(c, a, p);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    fxt_matll_permrow(a, p);

    a0 = fxt_matll_clone(a, 0, m-1, 0, m-1);
    a1 = fxt_matll_clone(a, m, n-1, 0, m-1);

    r = fxt_matll_new(n-m, m);
    cc = fxt_matld_to_matll(c);
    fxt_matll_mul(r, cc, a0);
    fxt_matll_sub(r, a1);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (i=0; i< n-m; i++)
      for (j=0; j< m; j++)
	if (fabs((double)MENT(r, i, j)) > 1e-13) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matll_ipstab: error in inversion\n");
	  return;
	}

    fxt_matll_del(cc);
    fxt_matll_del(r);
    fxt_matll_del(a1);
    fxt_matll_del(a0);
    fxt_matld_del(c);
    fxt_vecl_del(p);
    fxt_matll_del(a);

    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }
}
