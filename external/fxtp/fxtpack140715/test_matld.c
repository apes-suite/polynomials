#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"
#include "fxt_vecld.h"
#include "fxt_vecl.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing matld routines
*********************************************************************/

void test_matld(void) {
  fxt_matld *a, *b, *c;
  fxt_vecld *x, *y;
  long i, j, k, n, m, nm, n0, n1, m0, m1;
  double v, v0, v1, prec;

  printf("- test_matld...\n");

  for (n=0; n <= 10; n += 5)
    for (m=0; m <= 10; m += 5) {

      /* check create and delete */
      a = fxt_matld_new(n, m);

      for (i=0; i< n*m; i++)
	a->v[i] = 0.0;

      fxt_matld_del(a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      if (n > 0 && m > 0) {	

	/* check shrink */
	a = fxt_matld_new(n + 2, m + 2);

	for (n0 = 0; n0 <= 2; n0 += 2)
	  for (m0 = 0; m0 <= 2; m0 += 2) {
	    
	    for (i=0; i< n + 2; i++)
	      for (j=0; j< m + 2; j++)
		MENT(a, i, j) = (double) (i * 10 + j);

	    fxt_matld_shrink(a, n - n0, m - m0);
	    if (fxt_error_level() > FXT_ERROR_WARN)
	      return;

	    for (i=0; i< a->nrow; i++)
	      for (j=0; j < a->ncol; j++)
		if (MENT(a, i, j) != (double) (i * 10 + j)) {
		  fxt_error_set(FXT_ERROR_FXTBUG,
				"test_matld: error in shrink\n");
		  return;
		}

	    a->nrow = n + 2;
	    a->ncol = m + 2;
	  }

	fxt_matld_del(a);

	/* check clone */
	a = fxt_matld_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	b = fxt_matld_clone(a, 0, n-1, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(b, i, j) != (double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in clone\n");
	      return;
	    }

	fxt_matld_del(b);

	n0 = n / 3;
	n1 = 2 * n / 3;
	m0 = m / 3;
	m1 = 2 * m / 3;

	b = fxt_matld_clone(a, n0, n1, m0, m1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i= n0; i<= n1; i++)
	  for (j= m0; j<= m1; j++)
	    if (MENT(b, i-n0, j-m0) != (double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in clone\n");
	      return;
	    }

	/* check zeroing */
	fxt_matld_zero(b, 0, b->nrow-1, 0, b->ncol-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< b->nrow * b->ncol; i++)
	  if (b->v[i] != 0.0) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: error in zero\n");
	    return;
	  }

	fxt_matld_del(b);

	fxt_matld_zero(a, n0, n1, m0, m1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if ((i < n0 || n1 < i) || (j < m0 || m1 < j)) {
	      if (MENT(a, i, j) != (double) (i * 10 + j)) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in zero\n");
		return;
	      }
	    } else {
	      if (MENT(a, i, j) != 0.0) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in zero\n");
		return;
	      }
	    }

	fxt_matld_del(a);
      }	/* if (n > 0 && m > 0) */

      /* check unit matrix */
      if (n == m) {
	a = fxt_matld_new(n, n);

	fxt_matld_unit(a);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< n; j++)
	    if (MENT(a, i, j) != (i == j ? 1.0 : 0.0)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in unit\n");
	      return;
	    }

	fxt_matld_del(a);
      }	/* if n == m */

      /* check copy */
      if (n > 0 && m > 0) {
	a = fxt_matld_new(n + 10, m + 10);
	b = fxt_matld_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(b, i, j) = (double) (i * 10 + j);

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    MENT(a, i, j) = -1.0;

	fxt_matld_prtset(a, 5, 5+n-1, 5, 5+m-1, b);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    if (i < 5 || 5 + n <= i || j < 5 || 5 + m <= j) {
	      if (MENT(a, i, j) != -1.0) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in prtset\n");
		return;
	      }
	    } else {
	      if (MENT(a, i, j) != (double) ((i-5) * 10 + (j-5))) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in prtset\n");
		return;
	      }
	    }
      
	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(b, i, j) = -1.0;

	fxt_matld_setprt(b, a, 5, 5+n-1, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(b, i, j) != (double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in setprt\n");
	      return;
	    }

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    MENT(a, i, j) = -1.0;

	fxt_matld_prtsetprt(a, 6, 5+n-2, 6, 5+m-2, b, 1, n-2, 1, m-2);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n+10; i++)
	  for (j=0; j< m+10; j++)
	    if (i < 6 || 4 + n <= i || j < 6 || 4 + m <= j) {
	      if (MENT(a, i, j) != -1.0) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in prtsetprt\n");
		return;
	      }
	    } else {
	      if (MENT(a, i, j) != (double) ((i-5) * 10 + (j-5))) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in prtsetprt\n");
		return;
	      }
	    }

	fxt_matld_del(b);
	fxt_matld_del(a);
      }

      /* check subtraction */
      a = fxt_matld_new(n, m);
      b = fxt_matld_new(n, m);

      for (i=0; i< n; i++)
	for (j=0; j< m; j++) {
	  MENT(a, i, j) = (double) (i * 10 + j);
	  MENT(b, i, j) = (double) (i * 10 + j + 1);
	}

      fxt_matld_sub(b, a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n * m; i++)
	if (b->v[i] != 1.0) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld: error in sub\n");
	  return;
	}

      fxt_matld_del(b);
      fxt_matld_del(a);

      /* check matrix multiply */
      for (nm = 0; nm <= 10; nm += 5) {
	a = fxt_matld_new(n, nm);
	b = fxt_matld_new(nm, m);
	c = fxt_matld_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< nm; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	for (i=0; i< nm; i++)
	  for (j=0; j< m; j++)
	    MENT(b, i, j) = (double) (i * 10 + j);

	fxt_matld_mul(c, a, b);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++) {
	    v = 10.0 * (double) (nm * i * j);

	    for (k = 0; k < nm; k++)
	      v += 10.0 * k * k + (100 * i + j) * k;

	    if (MENT(c, i, j) != v) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in mul\n");
	      return;
	    }
	  }

	fxt_matld_del(c);
	fxt_matld_del(b);
	fxt_matld_del(a);
      }	/* for nm */

      /* Frobenius norm */
      a = fxt_matld_new(n, m);

      v = 0.0;
      for (i=0; i< n; i++)
	for (j=0; j< m; j++) {
	  MENT(a, i, j) = (double) (i * 10 + j);
	  v += (double) (i * 10 + j) * (double) (i * 10 + j);
	}

      if (v == 0.0)
	v = fxt_matld_normf(a);
      else
	v = fxt_matld_normf(a) / sqrt(v) - 1.0;
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      if (fabs(v) > 1e-15) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_matld: error in normf\n");
	return;
      }

      /* 2-norm */
      v = fxt_matld_normf(a);
      v0 = 0.0;
      for (prec = 1e-12; prec < 1e-10; prec *= 10.0) {
	v1 = fxt_matld_norm2(a, prec);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	/* 2norm <= Fnorm <= min(n,m) * 2norm */
	if (v1 * (1.0 - prec) > v ||
	    v > v1 * (1.0 + prec) * (n < m ? n : m)) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld: error in norm2\n");
	  return;
	}

	/* convergence to the same value */
	if (v0 != 0.0) {
	  if (v1 * (1.0 - prec) > v0 || v0 > v1 * (1.0 + prec)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: convergence error in norm2\n");
	    return;
	  }
	}

	v0 = v1;
      }

      fxt_matld_del(a);

      if (n > 0 && m > 0) {
	/* check swap */
	fxt_vecl *p;

	a = fxt_matld_new(n, m);

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	p = fxt_vecl_new(m);

	for (i=0; i< m; i++)
	  p->v[i] = i;

	for (nm = 0; nm < 20; nm ++) {
	  m0 = (long) (m * (double) rand() / (1.0 + RAND_MAX));
	  m1 = (long) (m * (double) rand() / (1.0 + RAND_MAX));

	  if (0 > m0 || m0 >= m || 0 > m1 || m1 >= m) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: irregal index?\n");
	    return;
	  }

	  fxt_vecl_swap(p, m0, m1);
	  fxt_matld_swapcol(a, m0, m1);
	  if (fxt_error_level() > FXT_ERROR_WARN)
	    return;

	  for (i=0; i< n; i++)
	    for (j=0; j< m; j++)
	      if (MENT(a, i, j) != (double) (i * 10 + p->v[j])) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in swapcol\n");
		return;
	      }
	}

	/* check permute columns */
	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	fxt_matld_permcol(a, p);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (double) (i * 10 + p->v[j])) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in permcol\n");
	      return;
	    }
	
	fxt_vecl_del(p);

	/* check permute rows */
	p = fxt_vecl_new(n);

	for (i=0; i< n; i++)
	  p->v[i] = i;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	for (nm=0; nm < 20; nm ++) {
	  n0 = (long) (n * (double) rand() / (1.0 * RAND_MAX));
	  n1 = (long) (n * (double) rand() / (1.0 * RAND_MAX));

	  if (0 > n0 || n0 >= n || 0 > n1 || n1 >= n) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: irregal index?\n");
	    return;
	  }

	  fxt_vecl_swap(p, n0, n1);
	  fxt_matld_swaprow(a, n0, n1);
	  if (fxt_error_level() > FXT_ERROR_WARN)
	    return;

	  for (i=0; i< n; i++)
	    for (j=0; j< m; j++)
	      if (MENT(a, i, j) != (double) p->v[i] * 10 + j) {
		fxt_error_set(FXT_ERROR_FXTBUG,
			      "test_matld: error in swaprow\n");
		return;
	      }
	}

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    MENT(a, i, j) = (double) (i * 10 + j);

	fxt_matld_permrow(a, p);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (double) (p->v[i] * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in permrow\n");
	      return;
	    }

	fxt_matld_ipermrow(a, p);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  for (j=0; j< m; j++)
	    if (MENT(a, i, j) != (double) (i * 10 + j)) {
	      fxt_error_set(FXT_ERROR_FXTBUG,
			    "test_matld: error in ipermrow\n");
	      return;
	    }

	fxt_vecl_del(p);
	fxt_matld_del(a);
      }	/* if (n > 0 && m > 0) */

      /* check get row/column */
      a = fxt_matld_new(n, m);
      x = fxt_vecld_new(m);
      y = fxt_vecld_new(n);

      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  MENT(a, i, j) = (double) (i * 10 + j);

      for (i=0; i< n; i++) {
	fxt_matld_getrow(x, a, i);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (i * 10 + j)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: error in getrow\n");
	    return;
	  }
      }

      for (j=0; j< m; j++) {
	fxt_matld_getcol(y, a, j);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 + j)) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: error in getcol\n");
	    return;
	  }
      }

      /* check 2-norms */

      fxt_matld_norm2rows(y, a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++) {
	v = 0.0;
	for (j=0; j< m; j++)
	  v += (double) (i * 10 + j) * (double) (i * 10 + j);

	if (v == 0.0)
	  v = y->v[i] * y->v[i];
	else
	  v = y->v[i] * y->v[i] / v - 1.0;

	if (fabs(v) > 1e-15) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld: error in norm2rows\n");
	  return;
	}
      }

      fxt_matld_norm2cols(x, a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (j=0; j< m; j++) {
	v = 0.0;
	for (i=0; i< n; i++)
	  v += (double) (i * 10 + j) * (double) (i * 10 + j);

	if (v == 0.0)
	  v = x->v[j] * x->v[j];
	else
	  v = x->v[j] * x->v[j] / v - 1.0;

	if (fabs(v) > 1e-15) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld: error in norm2cols\n");
	  return;
	}
      }

      /* check scale rows/cols */

      if (n > 0 && m > 0) {
	fxt_vecld_einv(x);
	fxt_matld_scalecol(a, x);
	fxt_matld_norm2cols(x, a);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (fabs(x->v[j] - 1.0) > 1e-15) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: error in scalecol\n");
	    return;
	  }

	fxt_matld_norm2rows(y, a);
	fxt_vecld_einv(y);
	fxt_matld_scalerow(a, y);
	fxt_matld_norm2rows(y, a);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (fabs(y->v[i] - 1.0) > 1e-15) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld: error in scalerow\n");
	    return;
	  }
      } else {
	/* just call and check size */
	fxt_matld_scalecol(a, x);
	fxt_matld_scalerow(a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
      }

      fxt_vecld_del(y);
      fxt_vecld_del(x);
      fxt_matld_del(a);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }
}
