#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test vector of double
*********************************************************************/

static int differ(double *x, double *y, long length) {
  long i;

  for (i=0; i< length; i++)
    if (x[i] != y[i])
      return 1;

  return 0;
}

void test_vecld(void) {
  fxt_vecld *x, *y;
  double s;			/* norm */
  fxt_vecl *p;			/* perm */
  long i, n;

  printf("- test_vecld...\n");

  /* create and delete */
  for (n=0; n <= 10; n += 10) {
    x = fxt_vecld_new(n);

    if (x->n != n) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in new\n");
      return;
    }

    for (i=0; i< n; i++)
      x->v[i] = 0.0;

    fxt_vecld_del(x);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  /* clone */
  x = fxt_vecld_new(10);

  for (i=0; i< 10; i++)
    x->v[i] = (double) i;

  y = fxt_vecld_clone(x, 0, 9);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 10)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecld: error in clone\n");
    return;
  }

  fxt_vecld_del(y);

  y = fxt_vecld_clone(x, 2, 7);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v + 2, y->v, 6)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecld: error in clon\n");
    return;
  }

  fxt_vecld_del(y);

  /* check all zero */
  fxt_vecld_zero(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != 0.0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in zero\n");
      return;
    }

  /* check all one */
  fxt_vecld_one(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != 1.0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in one\n");
      return;
    }

  /* check random */
  fxt_vecld_rand(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (0 > x->v[i] || x->v[i] >= 1.0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in rand\n");
    }

  /* check unit vector */
  for (n=0; n< 10; n++) {
    fxt_vecld_unit(x, n);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (i=0; i< 10; i++)
      if ((i == n && x->v[i] != 1.0) ||
	  (i != n && x->v[i] != 0.0)) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_vecld: error in unit\n");
	return;
      }
  }

  /* check negation */
  for (i=0; i< 10; i++)
    x->v[i] = (double) i;

  fxt_vecld_neg(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != (double) -i) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in neg\n");
      return;
    }

  /* check inversion */
  for (i=0; i< 10; i++)
    x->v[i] = (double) (i + 1);

  fxt_vecld_einv(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    double d;

    d = x->v[i] * (i + 1) - 1.0;

    if (fabs(d) > 2e-16 * (i+1)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in einv\n");
      return;
    }
  }

  /* check square root */
  for (i=0; i< 10; i++)
    x->v[i] = (double) i;

  fxt_vecld_sqrt(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    double d;

    d = x->v[i] * x->v[i] - (double) i;

    if (i != 0)
      d /= i;

    if (fabs(d) > 1e-15) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in sqrt\n");
      return;
    }
  }

  /* check square */
  fxt_vecld_sq(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    double d;

    d = x->v[i] - (double) i;
    if (i != 0)
      d /= i;

    if (fabs(d) > 1e-15) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in sq\n");
      return;
    }
  }

  /* check scale */
  for (i=0; i< 10; i++)
    x->v[i] = (double) (i + 1);

  fxt_vecld_scale(x, 2, 7, 3.0);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    double d, z;

    z = ((2 <= i && i <= 7) ? 3.0 : 1.0) * (i + 1);

    d = (x->v[i] - z) / z;

    if (fabs(d) > 1e-15) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in scale\n");
      return;
    }
  }

  /* check add */
  for (i=0; i< 10; i++)
    x->v[i] = (double) i;

  y = fxt_vecld_clone(x, 0, 9);

  fxt_vecld_add(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != (double) (i * 2)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in add\n");
      return;
    }

  /* check subtract */
  fxt_vecld_sub(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 10)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecld: error in sub\n");
    return;
  }

  /* check multiply */
  fxt_vecld_emul(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != (double) (i * i)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in emul\n");
      return;
    }

  /* check dot product */
  s = fxt_vecld_norm2(y);

  s *= s;
  for (i=0; i< 10; i++)
    s -= (double) (i * i);

  if (fabs(s) > 1e-13) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecld: error in norm2\n");
    return;
  }

  /* check normalize */
  fxt_vecld_normalize(x);

  if (fabs(1.0 - fxt_vecld_norm2(x)) > 1e-15) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecld: error in normalize\n");
    return;
  }

  /* check permutation */
  p = fxt_vecl_new(10);

  for (i=0; i< 10; i++)
    p->v[i] = i;

  for (n=0; n< 20; n++) {
    long j, k;

    j = (long) (10.0 * rand() / (1.0 + RAND_MAX));
    k = (long) (10.0 * rand() / (1.0 + RAND_MAX));

    if (0 > j || j >= 10 || 0 > k || k >= 10) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: irregal index?\n");
      return;
    }

    fxt_vecl_swap(p, j, k);
  }

  fxt_vecld_perm(y, p);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (y->v[i] != (double) p->v[i]) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in perm\n");
      return;
    }

  /* check inverse permutation */
  fxt_vecld_iperm(y, p);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (y->v[i] != (double) i) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecld: error in iperm\n");
      return;
    }

  fxt_vecl_del(p);
  fxt_vecld_del(y);
  fxt_vecld_del(x);
}
