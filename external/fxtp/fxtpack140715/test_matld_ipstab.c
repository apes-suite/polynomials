#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test stable interpolation for matld
*********************************************************************/

void test_matld_ipstab(void) {
  long n, m, i, j;
  fxt_matld *a, *c, *a0, *a1, *r;
  fxt_vecl *p;
  int cc;

  printf("- test_matld_ipstab...\n");

  for (cc = 0; cc < 2; cc ++) {
    if (cc == 0) {
      n = 50;  m = 20;
    } else {
      n = 50;  m = 30;
    }

    a = fxt_matld_new(n, m);
    for (i=0; i< n*m; i++)
      a->v[i] = (double) (rand() / (1.0 + RAND_MAX));

    p = fxt_vecl_new(n);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    c = fxt_matld_ipstab(a, p);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    fxt_matld_permrow(a, p);

    a0 = fxt_matld_clone(a, 0, m-1, 0, m-1);
    a1 = fxt_matld_clone(a, m, n-1, 0, m-1);

    r = fxt_matld_new(n-m, m);
    fxt_matld_mul(r, c, a0);
    fxt_matld_sub(r, a1);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (i=0; i< n-m; i++)
      for (j=0; j< m; j++)
	if (fabs(MENT(r, i, j)) > 1e-13) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_ipstab: inversion error\n");
	  return;
	}

    fxt_matld_del(r);
    fxt_matld_del(a1);
    fxt_matld_del(a0);
    fxt_matld_del(c);
    fxt_vecl_del(p);
    fxt_matld_del(a);

    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }
}
