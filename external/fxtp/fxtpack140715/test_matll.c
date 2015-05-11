#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matll.h"
#include "fxt_vecll.h"
#include "fxt_vecl.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test for matrix of long double
*********************************************************************/

void test_matll(void) {
  fxt_matll *a, *b, *c;
  fxt_vecll *x, *y;
  long i, j, k, n, m, nm, n0, n1, m0, m1;
  long double v;
  fxt_matld *d;

  printf("- test_matll...\n");

  for (n=0; n <= 10; n += 5)
    for (m=0; m <= 10; m += 5) {

      /* check create and delete */
      a = fxt_matll_new(n, m);

      for (i=0; i< n*m; i++)
	a->v[i] = 0.0;

      fxt_matll_del(a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      if (n > 0 && m > 0) {	

	/* check clone */
	a = fxt_matll_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (long double) (i * 10 + j);

	b = fxt_matll_clone(a, 0, n-1, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(b, i, j) != (long double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in clone\n");
	      return;
	    }

	fxt_matll_del(b);

	n0 = n / 3;
	n1 = 2 * n / 3;
	m0 = m / 3;
	m1 = 2 * m / 3;

	b = fxt_matll_clone(a, n0, n1, m0, m1);

	for (i= n0; i<= n1; i++)
	  for (j= m0; j<= m1; j++)
	    if (MENT(b, i-n0, j-m0) != (long double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in clone\n");
	      return;
	    }

	fxt_matll_del(b);
	fxt_matll_del(a);
      }	/* if (n > 0 && m > 0) */

      /* check copy */
      if (n > 0 && m > 0) {
	a = fxt_matll_new(n + 10, m + 10);
	b = fxt_matll_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(b, i, j) = (long double) (i * 10 + j);

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    MENT(a, i, j) = -1.0;

	fxt_matll_set(a, 5, 5+n-1, 5, 5+m-1, b);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    if (i < 5 || 5 + n <= i || j < 5 || 5 + m <= j) {
	      if (MENT(a, i, j) != -1.0) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matll: error in set\n");
		return;
	      }
	    } else {
	      if (MENT(a, i, j) !=
		  (long double) ((i-5) * 10 + (j-5))) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matll: error in set\n");
		return;
	      }
	    }

	fxt_matll_del(b);
	fxt_matll_del(a);
      }

      /* check subtraction */
      a = fxt_matll_new(n, m);
      b = fxt_matll_new(n, m);

      for (i=0; i< n; i++)
	for (j=0; j< m; j++) {
	  MENT(a, i, j) = (long double) (i * 10 + j);
	  MENT(b, i, j) = (long double) (i * 10 + j + 1);
	}

      fxt_matll_sub(b, a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n * m; i++)
	if (b->v[i] != 1.0) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matll: error in sub\n");
	  return;
	}
      
      fxt_matll_del(b);
      fxt_matll_del(a);

      /* check matrix multiply */
      for (nm = 0; nm <= 10; nm += 5) {
	a = fxt_matll_new(n, nm);
	b = fxt_matll_new(nm, m);
	c = fxt_matll_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< nm; j++)
	    MENT(a, i, j) = (long double) (i * 10 + j);

	for (i=0; i< nm; i++)
	  for (j=0; j< m; j++)
	    MENT(b, i, j) = (long double) (i * 10 + j);

	fxt_matll_mul(c, a, b);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++) {
	    v = 10.0 * (long double) (nm * i * j);

	    for (k = 0; k < nm; k++)
	      v += 10.0 * k * k + (100 * i + j) * k;

	    if (MENT(c, i, j) != v) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in mul\n");
	      return;
	    }
	  }

	fxt_matll_del(c);
	fxt_matll_del(b);
	fxt_matll_del(a);
      }	/* for nm */

      /* Frobenius norm */
      a = fxt_matll_new(n, m);

      v = 0.0;
      for (i=0; i< n; i++)
	for (j=0; j< m; j++) {
	  MENT(a, i, j) = (long double) (i * 10 + j);
	  v += (long double) (i*10 + j) * (long double) (i*10 + j);
	}

      if (v == 0.0)
	v = fxt_matll_normf(a);
      else
	v = fxt_matll_normf(a) / sqrt(v) - 1.0;
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      if (fabs(v) > 1e-15) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_matll: error in normf\n");
	return;
      }

      fxt_matll_del(a);

      if (n > 0 && m > 0) {
	/* check permute rows */
	fxt_vecl *p;

	a = fxt_matll_new(n, m);

	p = fxt_vecl_new(n);

	for (i=0; i< n; i++)
	  p->v[i] = i;

	for (nm=0; nm < 20; nm ++) {
	  n0 = (long) (n * (long double) rand() / (1.0 * RAND_MAX));
	  n1 = (long) (n * (long double) rand() / (1.0 * RAND_MAX));

	  if (0 > n0 || n0 >= n || 0 > n1 || n1 >= n) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matll: irregal index?\n");
	    return;
	  }

	  fxt_vecl_swap(p, n0, n1);
	}

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (long double) (i * 10 + j);

	fxt_matll_permrow(a, p);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (long double) (p->v[i] * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in permrow\n");
	      return;
	    }

	fxt_matll_ipermrow(a, p);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (long double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in ipermrow\n");
	      return;
	    }

	fxt_vecl_del(p);
	fxt_matll_del(a);
      }	/* if (n > 0 && m > 0) */

      /* check get row/column */
      a = fxt_matll_new(n, m);
      x = fxt_vecll_new(m);
      y = fxt_vecll_new(n);

      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  MENT(a, i, j) = (long double) (i * 10 + j);

      for (i=0; i< n; i++) {
	fxt_matll_getrow(x, a, i);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (long double) (i * 10 + j)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matll: error in getrow\n");
	    return;
	  }
      }

      for (j=0; j< m; j++) {
	fxt_matll_getcol(y, a, j);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (long double) (i * 10 + j)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matll: error in getcol\n");
	    return;
	  }
      }

      /* check scale rows/cols */

      if (n > 0 && m > 0) {
	for (i=0; i< n; i++)
	  y->v[i] = (double) i;

	fxt_matll_scalerow(a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (long double) ((i * 10 + j) * i)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matll: error in scalerow\n");
	      return;
	    }
      } else {
	/* just call and check size */
	fxt_matll_scalerow(a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
      }

      fxt_vecll_del(y);
      fxt_vecll_del(x);

      /* check conversions */
      d = fxt_matll_to_matld(a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  if (MENT(d, i, j) != (double) ((i * 10 + j) * i)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matll: error in matll_to_matld\n");
	    return;
	  }

      b = fxt_matld_to_matll(d);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  if (MENT(b, i, j) != (long double) ((i * 10 + j) * i)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matll: error in matld_to_matll\n");
	    return;
	  }

      fxt_matll_del(b);
      fxt_matld_del(d);

      fxt_matll_del(a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }
}
