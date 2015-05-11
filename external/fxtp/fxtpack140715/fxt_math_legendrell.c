#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_math.h"

#include "fxt_vecll.h"
#include "fxt_matll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  associated Legendre functions in long double
*********************************************************************/


/*** extended exponent ***/
typedef struct {
  long e;			/* exponent */
  long double m;		/* mantissa */
} xx;

/*** initialization ***/
static void xx_init(void);

/*** set value ***/
static void xx_setval(xx *xq, long double v);

/*** get value ***/
static long double xx_getval(xx *xq);

/*** multiplication: x1 *= v2 ***/
static void xx_mull(xx *x1, long double v2);

/*** x1 = x1 * v2 + x3 * v4 ***/
static void xx_rec(xx *x1, long double v2, xx *x3, long double v4);



/*** func(0, 0, x) ***/
static long double func00(long double x) {
  return 1.0L;
}

static long double coef0(long m) {
  return fxt_lsqrt((long double)(2 * m + 1) / (long double)(2 * m));
}

static long double coef1(long m) {
  return fxt_lsqrt((long double)(2 * m + 3));
}

static long double coef2(long m, long n) {
  long nn = (n + 1) * (n + 1);
  return fxt_lsqrt((long double)(4*nn-1) / (long double)(nn - m*m));
}

static long double coef3(long m, long n) {
  return fxt_lsqrt((long double)(n*n-m*m) / (long double)(4*n*n-1));
}


/*** compute the function directly ***/
static xx funcmn(long m, long n, long double x) {
  xx p, p0, p1;
  long double qx;
  long i;

  /* qx = sin(th) */
  qx = fxt_lsqrt(1.0L - x * x);

  /* p = func(0, 0, x) */
  xx_setval(&p, func00(x));

  /* p = func(m, m, x) */
  for (i=1; i<= m; i++)
    xx_mull(&p, coef0(i) * qx);

  if (n == m)
    return p;

  /* p0 = func(m, m, x) */
  p0 = p;

  /* p1 = func(m, m+1, x) */
  p1 = p;
  xx_mull(&p1, coef1(m) * x);

  /* p1 = func(m, i+1, x) */
  for (i= m+1; i< n; i++) {
    p = p0;
    p0 = p1;

    /* recursion */
    xx_rec(&p1, coef2(m, i) * x, &p, -coef2(m, i) * coef3(m, i));
  }

  /* p1 = func(m, n, x) */
  return p1;
}


/*** compute function matrix ***/
fxt_matll* fxt_legendre_matll(fxt_vecll *gv, long m, long nmax,
			      long odd) {
  fxt_matll *mat;		/* function matrix */
  long nx;			/* number of points */
  long nn;			/* number of function orders */
  long i, j, n;
  
  /* check null pointers */
  if (gv == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_legendre_matll: null pointer\n");
    return NULL;
  }

  /* check input */
  if (m + odd > nmax) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_legendre_matll: irregal index\n");
    return NULL;
  }

  /* initialize extended exponent */
  xx_init();

  /* number of points */
  nx = (gv->n + 1 - odd) / 2;

  /* matrix size */
  nn = (nmax - m + 1 + 1 - odd) / 2;

  /* allocate matrix */
  mat = fxt_matll_new(nx, nn);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* compute for all points... */
  for (i=0; i< nx; i++) {
    xx p0, p1, p;

    /* the point */
    long double x = gv->v[i];

    /* squared x */
    long double x2 = x * x;

    /* the first value */
    p0 = funcmn(m, m + odd, x);
    MENT(mat, i, 0) = xx_getval(&p0);

    /* the second value */
    if (nn > 1) {
      p1 = funcmn(m, m + 2 + odd, x);
      MENT(mat, i, 1) = xx_getval(&p1);
    }

    /* the other values */
    for (j=2, n= m+odd; n+4 <= nmax; j++, n+=2) {
      long double a = coef2(m, n+2) * coef2(m, n+3);

      long double b = -coef2(m, n+3) * coef3(m, n+3)
	- a * coef3(m, n+2) / coef2(m, n+1);

      long double c = - a * coef3(m, n+2) * coef3(m, n+1);

      p = p0;
      p0 = p1;
      
      xx_rec(&p1, a * x2 + b, &p, c);
      MENT(mat, i, j) = xx_getval(&p1);
    }
  }

  return mat;
}

/*********************************************************************
  extended exponent
*********************************************************************/

/*** constants for extended exponent ***/
static long double qbase, qibase, qhbase, qihbase;


/*** initialization ***/
static void xx_init(void) {
  static int exec = 0;

  /* check initialized */
  if ((exec++) != 0)
    return;

  /* qhbase = 2^128 */
  qhbase = (long double) ldexp(1.0, 128);

  /* qbase = 2^256 */
  qbase = qhbase * qhbase;

  /* qibase = 2^{-256} */
  qibase = 1.0L / qbase;

  /* qihbase = 2^{-128} */
  qihbase = 1.0L / qhbase;
}


/*** force the exponent ***/
static void xx_setexp(xx *xq, long ex) {
  /* scale down */
  while (xq->e > ex) {
    xq->e --;
    xq->m *= qbase;
  }

  /* scale up */
  while (xq->e < ex) {
    xq->e ++;
    xq->m *= qibase;
  }
}


/*** commonize the exponent ***/
static void xx_comexp(xx *x1, xx *x2) {
  /* check zeros */
  if (x1->m == 0.0 || x2->m == 0.0)
    return;

  /* commonize exponent */
  if (x1->e > x2->e)
    xx_setexp(x2, x1->e);
  else if (x2->e > x1->e)
    xx_setexp(x1, x2->e);
}


/*** normalize ***/
static void xx_normal(xx *xq) {
  int negative;			/* sign bit */
  long double m;		/* absolute mantissa */

  /* check zero */
  if (xq->m == 0.0) {
    xq->e = 0;
    return;
  }

  /* remove sign bit */
  negative = (xq->m < 0.0);
  m = (negative ? -(xq->m) : xq->m);

  /* normalize */
  while (m >= qhbase) {
    xq->e ++;
    m *= qibase;
  }

  while (m < qihbase) {
    xq->e --;
    m *= qbase;
  }

  /* restore sign bit */
  xq->m = negative ? -m : m;
}


/*** set value ***/
static void xx_setval(xx *xq, long double v) {
  /* set value */
  xq->e = 0;
  xq->m = v;

  /* normalize */
  xx_normal(xq);
}


/*** get value ***/
static long double xx_getval(xx *xq) {
  /* make a copy */
  xx x = *xq;

  /* force exponent to be zero */
  xx_setexp(&x, 0);

  /* return value */
  return x.m;
}


/*** multiplication: x1 *= v2 ***/
static void xx_mull(xx *x1, long double v2) {
  /* multiply */
  x1->m *= v2;

  /* normalize */
  xx_normal(x1);
}


/*** x1 = x1 * v2 + x3 * v4 ***/
static void xx_rec(xx *x1, long double v2, xx *x3, long double v4) {
  xx c;

  /* c = x1 * x2 */
  x1->m *= v2;
  xx_normal(x1);

  /* check for zero */
  if (x3->m == 0.0 || v4 == 0.0)
    return;

  /* c = x3 * v4 */
  c = *x3;
  c.m *= v4;
  xx_normal(&c);

  /* add c to x1 */
  xx_comexp(x1, &c);
  x1->m += c.m;
  xx_normal(x1);
}
