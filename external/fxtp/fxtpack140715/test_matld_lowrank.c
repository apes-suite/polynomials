#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test lowrank approximate decomposition for matld
*********************************************************************/

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void test_matld_lowrank(void) {
  long n, m, i, j, k, prev_k;
  fxt_vecld *x, *y;
  fxt_matld *a, *(z[2]), *aa;
  int errc;

  printf("- test_matld_lowrank...\n");

  for (n = 25; n <= 50; n += 25)
    for (m = 25; m <= 50; m += 25) {
      double eps = 0.5;		/* precision */

      x = fxt_vecld_new(n);
      y = fxt_vecld_new(m);

      fxt_vecld_rand(x);	/* [0, 1] */
      fxt_vecld_rand(y);	/* [2, 3] */
      for (i=0; i< y->n; i++)
	y->v[i] += 2.0;
      
      a = fxt_matld_new(n, m);
      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  MENT(a, i, j) = 1.0 / (x->v[i] - y->v[j]);
      
      prev_k = 0;
      do {
	z[0] = fxt_matld_new(n, MIN(n, m));
	z[1] = fxt_matld_new(MIN(n, m), m);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	errc = fxt_matld_lowrank(a, eps, z[0], z[1]);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	if (errc != 0) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_lowrank: rank overflow\n");
	  return;
	}

	k = z[0]->ncol;

	aa = fxt_matld_new(n, m);
	fxt_matld_mul(aa, z[0], z[1]);
	fxt_matld_sub(aa, a);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	if (k < prev_k) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_lowrank: abnormal rank\n");
	  return;
	}

	if (fxt_matld_normf(aa)
	    > eps + MAX(a->nrow, a->ncol) * 2e-16) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_lowrank: decomposition error\n");
	  return;
	}

	fxt_matld_del(aa);
	fxt_matld_del(z[1]);
	fxt_matld_del(z[0]);

	eps *= 0.1;
	prev_k = k;

	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
      } while (eps > 1e-15 && k < MIN(n, m));

      fxt_matld_del(a);
      fxt_vecld_del(y);
      fxt_vecld_del(x);

      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }
}
