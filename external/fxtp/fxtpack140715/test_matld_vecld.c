#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test matrix-vector product for matld
*********************************************************************/

static long sum(long n0, long n1) {
  long n, s;

  s = 0;

  for (n = n0; n < n1; n ++)
    s += n;

  return s;
}

static long sumsq(long n0, long n1) {
  long n, s;

  s = 0;

  for (n = n0; n < n1; n ++)
    s += n * n;

  return s;
}

void test_matld_vecld(void) {
  long i, j, n, m;
  fxt_matld *a;
  fxt_vecld *x, *y, *xx, *yy;
  
  printf("- test_matld_vecld...\n");

  for (n=0; n<= 10; n+=5)
    for (m=0; m<= 10; m+=5) {
      a = fxt_matld_new(n, m);
      x = fxt_vecld_new(m);
      y = fxt_vecld_new(n);
      xx = fxt_vecld_new(m + 10);
      yy = fxt_vecld_new(n + 10);

      /* check setmatvec */
      for (i=0; i< n; i++)
	for (j=0; j< m; j++)
	  MENT(a, i, j) = (double) (i * 10 + j);

      for (j=0; j< m+10; j++)
	xx->v[j] = -1.0;

      for (j=0; j< m; j++)
	xx->v[j + 5] = x->v[j] = (double) j;

      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      fxt_matld_setmatvec(y, a, x);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	if (y->v[i] != (double) (i * 10 * sum(0, m)
				 + sumsq(0, m))) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_vecld: error in setmatvec\n");
	  return;
      }

      fxt_matld_addmatvec(y, a, x);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (i=0; i< n; i++)
	if (y->v[i] != (double) (i * 20 * sum(0, m)
				    + 2 * sumsq(0, m))) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_vecld: error in addmatvec\n");
	  return;
      }

      if (m > 0) {
	fxt_matld_setmatprtvec(y, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 * sum(0, m)
				      + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in setmatprtvec\n");
	    return;
	}

	fxt_matld_addmatprtvec(y, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
	
	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 20 * sum(0, m)
				      + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in addmatprtvec\n");
	    return;
	}

	fxt_matld_setmatprtvec(y, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 * sum(0, m)
				      + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in setmatprtvec\n");
	    return;
	}

	fxt_matld_addmatprtvec(y, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
	
	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 20 * sum(0, m)
				      + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in addmatprtvec\n");
	    return;
	}
      }
      
      if (n > 0) {
	fxt_matld_prtsetmatvec(y, 0, n-1, a, x);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 * sum(0, m)
				      + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatvec\n");
	    return;
	}

	fxt_matld_prtaddmatvec(y, 0, n-1, a, x);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 20 * sum(0, m)
				      + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatvec\n");
	    return;
	}

	fxt_matld_prtsetmatvec(yy, 5, 5+n-1, a, x);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 10 * sum(0, m)
					 + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatvec\n");
	    return;
	}

	fxt_matld_prtaddmatvec(yy, 5, 5+n-1, a, x);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 20 * sum(0, m)
					 + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatvec\n");
	    return;
	}
      }

      if (n > 0 && m > 0) {
	fxt_matld_prtsetmatprtvec(y, 0, n-1, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 * sum(0, m)
				      + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatprtvec\n");
	    return;
	}

	fxt_matld_prtaddmatprtvec(y, 0, n-1, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 20 * sum(0, m)
				      + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatprtvec\n");
	    return;
	}

	fxt_matld_prtsetmatprtvec(yy, 5, 5+n-1, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 10 * sum(0, m)
					 + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatprtvec\n");
	    return;
	}

	fxt_matld_prtaddmatprtvec(yy, 5, 5+n-1, a, x, 0, m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 20 * sum(0, m)
					 + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatprtvec\n");
	    return;
	}

	fxt_matld_prtsetmatprtvec(y, 0, n-1, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 10 * sum(0, m)
				      + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatprtvec\n");
	    return;
	}

	fxt_matld_prtaddmatprtvec(y, 0, n-1, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (y->v[i] != (double) (i * 20 * sum(0, m)
				      + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatprtvec\n");
	    return;
	}

	fxt_matld_prtsetmatprtvec(yy, 5, 5+n-1, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 10 * sum(0, m)
					 + sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsetmatprtvec\n");
	    return;
	}

	fxt_matld_prtaddmatprtvec(yy, 5, 5+n-1, a, xx, 5, 5+m-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (i=0; i< n; i++)
	  if (yy->v[i+5] != (double) (i * 20 * sum(0, m)
					 + 2 * sumsq(0, m))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddmatprtvec\n");
	    return;
	}
      }

      for (i=0; i< n+10; i++)
	yy->v[i] = -1;

      for (i=0; i< n; i++)
	yy->v[i+5] = y->v[i] = (double) i;

      fxt_matld_settrmatvec(x, a, y);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (j=0; j< m; j++)
	if (x->v[j] != (double) (10 * sumsq(0, n)
				    + j * sum(0, n))) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_vecld: "
			"error in settrmatvec\n");
	  return;
      }

      fxt_matld_addtrmatvec(x, a, y);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;

      for (j=0; j< m; j++)
	if (x->v[j] != (double) (20 * sumsq(0, n)
				    + 2 * j * sum(0, n))) {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_matld_vecld: "
			"error in addtrmatvec\n");
	  return;
      }

      if (n > 0) {
	fxt_matld_settrmatprtvec(x, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (10 * sumsq(0, n)
				      + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in settrmatprtvec\n");
	    return;
	}

	fxt_matld_addtrmatprtvec(x, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (20 * sumsq(0, n)
				      + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in addtrmatprtvec\n");
	    return;
	}

	fxt_matld_settrmatprtvec(x, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (10 * sumsq(0, n)
				      + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in settrmatprtvec\n");
	    return;
	}

	fxt_matld_addtrmatprtvec(x, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (20 * sumsq(0, n)
				      + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in addtrmatprtvec\n");
	    return;
	}
      }

      if (m > 0) {
	fxt_matld_prtsettrmatvec(x, 0, m-1, a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (10 * sumsq(0, n)
				      + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatvec\n");
	    return;
	}

	fxt_matld_prtaddtrmatvec(x, 0, m-1, a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (20 * sumsq(0, n)
				      + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatvec\n");
	    return;
	}

	fxt_matld_prtsettrmatvec(xx, 5, 5+m-1, a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (10 * sumsq(0, n)
					 + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatvec\n");
	    return;
	}

	fxt_matld_prtaddtrmatvec(xx, 5, 5+m-1, a, y);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (20 * sumsq(0, n)
					 + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatvec\n");
	    return;
	}	
      }

      if (n > 0 && m > 0) {
	fxt_matld_prtsettrmatprtvec(x, 0, m-1, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (10 * sumsq(0, n)
				      + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatprtvec\n");
	    return;
	}

	fxt_matld_prtaddtrmatprtvec(x, 0, m-1, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (20 * sumsq(0, n)
				      + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatprtvec\n");
	}

	fxt_matld_prtsettrmatprtvec(x, 0, m-1, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (10 * sumsq(0, n)
				      + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatprtvec\n");
	}

	fxt_matld_prtaddtrmatprtvec(x, 0, m-1, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (x->v[j] != (double) (20 * sumsq(0, n)
				      + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatprtvec\n");
	}

	fxt_matld_prtsettrmatprtvec(xx, 5, 5+m-1, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (10 * sumsq(0, n)
					 + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatprtvec\n");
	}

	fxt_matld_prtaddtrmatprtvec(xx, 5, 5+m-1, a, y, 0, n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (20 * sumsq(0, n)
					 + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatprtvec\n");
	}

	fxt_matld_prtsettrmatprtvec(xx, 5, 5+m-1, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (10 * sumsq(0, n)
					 + j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtsettrmatprtvec\n");
	}

	fxt_matld_prtaddtrmatprtvec(xx, 5, 5+m-1, a, yy, 5, 5+n-1);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;

	for (j=0; j< m; j++)
	  if (xx->v[j+5] != (double) (20 * sumsq(0, n)
					 + 2 * j * sum(0, n))) {
	    fxt_error_set(FXT_ERROR_FXTBUG,
			  "test_matld_vecld: "
			  "error in prtaddtrmatprtvec\n");
	}
      }

      fxt_vecld_del(yy);
      fxt_vecld_del(xx);
      fxt_vecld_del(y);
      fxt_vecld_del(x);
      fxt_matld_del(a);

      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }
}
