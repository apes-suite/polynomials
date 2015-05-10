#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_vecll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test vector of long double
*********************************************************************/

static int differ(long double *x, long double *y, long length) {
  long i;

  for (i=0; i< length; i++)
    if (x[i] != y[i])
      return 1;

  return 0;
}

static long double qabs(long double x) {
  if (x < 0)
    return -x;

  return x;
}

void test_vecll(void) {
  fxt_vecll *x, *y;
  long double eps;		/* machine epsilon */
  long double s;		/* norm */
  fxt_vecl *p;			/* perm */
  fxt_vecld *z;			/* conversion */
  long i, n;

  eps = 1.0L;

  while ((1.0L + eps) > 1.0L)
    eps *= 0.5L;

  printf("- test_vecll...\n");

  /* create and delete */
  for (n=0; n <= 10; n += 10) {
    x = fxt_vecll_new(n);

    if (x->n != n) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in new\n");
      return;
    }

    for (i=0; i< n; i++)
      x->v[i] = 0.0L;

    fxt_vecll_del(x);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  /* clone */
  x = fxt_vecll_new(10);

  for (i=0; i< 10; i++)
    x->v[i] = (long double) i;

  y = fxt_vecll_clone(x, 0, 9);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 10)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecll: error in clone\n");
    return;
  }

  fxt_vecll_del(y);

  y = fxt_vecll_clone(x, 2, 7);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v + 2, y->v, 6)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecll: error in clone\n");
    return;
  }

  fxt_vecll_del(y);

  /* check all one */
  fxt_vecll_one(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != 1.0L) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in one\n");
      return;
    }

  /* check negation */
  for (i=0; i< 10; i++)
    x->v[i] = (long double) i;

  fxt_vecll_neg(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != (long double) -i) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in neg\n");
      return;
    }

  /* check square root */
  for (i=0; i< 10; i++)
    x->v[i] = (long double) i;

  fxt_vecll_sqrt(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    long double d;

    d = x->v[i] * x->v[i] - (long double) i;

    if (i != 0)
      d /= i;

    if (qabs(d) > 4.0L * eps) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in sqrt\n");
      return;
    }
  }

  /* check square */
  fxt_vecll_sq(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    long double d;

    d = x->v[i] - (long double) i;
    if (i != 0)
      d /= i;

    if (qabs(d) > 4.0L * eps) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in sq\n");
      return;
    }
  }

  /* check scale */
  for (i=0; i< 10; i++)
    x->v[i] = (long double) (i + 1);

  fxt_vecll_scale(x, 2, 7, 3.0L);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++) {
    long double d, zz;

    zz = ((2 <= i && i <= 7) ? 3.0L : 1.0L) * (i + 1);

    d = (x->v[i] - zz) / zz;
    
    if (qabs(d) > 4.0L * eps) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in scale\n");
      return;
    }
  }

  /* check multiply*/
  for (i=0; i< 10; i++)
    x->v[i] = (long double) i;

  y = fxt_vecll_clone(x, 0, 9);

  fxt_vecll_emul(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (x->v[i] != (long double) (i * i)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in emul\n");
      return;
    }

  /* check dot product */
  s = fxt_vecll_dot(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    s -= (long double) (i * i * i);

  if (qabs(s) > 4.0L * eps) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecll: error in dot\n");
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
		    "test_vecll: irregal index?\n");
      return;
    }

    fxt_vecl_swap(p, j, k);
  }

  fxt_vecll_perm(y, p);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (y->v[i] != (long double) p->v[i]) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in perm\n");
      return;
    }

  /* check inverse permutation */
  fxt_vecll_iperm(y, p);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (y->v[i] != (long double) i) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in iperm\n");
      return;
    }

  fxt_vecl_del(p);
  fxt_vecll_del(y);

  /* check conversion to vecld */
  z = fxt_vecll_to_vecld(x);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  for (i=0; i< 10; i++)
    if (z->v[i] != (double) (i * i)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecll: error in vecll_to_vecld\n");
      return;
    }

  /* check conversion from vecld */
  y = fxt_vecld_to_vecll(z);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 10)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecll: error in vecld_to_vecll\n");
    return;
  }

  fxt_vecll_del(y);
  fxt_vecld_del(z);
  fxt_vecll_del(x);
}
