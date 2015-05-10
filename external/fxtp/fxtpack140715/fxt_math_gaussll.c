#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_math.h"

#include "fxt_vecll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  compute the Gauss nodes and weights
*********************************************************************/

static double pi = 3.14159265358979323846;

/*** compute Gauss nodes and weights ***/
void fxt_gauss_vecll(long n, fxt_vecll **xp, fxt_vecll **wp) {
  fxt_vecll *xv, *wv;		/* Gauss nodes and weights */
  long double eps;		/* machine epsilon */
  long nhalf;			/* half of n */
  long i, j, k;

  /* clear results for error return */
  *xp = *wp = NULL;

  /* check size */
  if (n <= 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_gauss_vecll: non-positive size\n");
    return;
  }

  /* machine epsilon */
  eps = 0.5L;
  while (1.0L + eps > 1.0L)
    eps *= 0.5L;
  eps *= 16.0L;

  /* a half of the number of the points */
  nhalf = (n + 1) / 2;

  /* the Gauss nodes */
  xv = fxt_vecll_new(n);
  if (fxt_error_raise() > FXT_ERROR_WARN || xv == NULL)
    return;

  /* the Gauss weights */
  wv = fxt_vecll_new(n);
  if (fxt_error_raise() > FXT_ERROR_WARN || xv == NULL)
    return;

  /* compute Gauss nodes: Pn(x) = 0 */
  for (i=0; i< nhalf; i++) {
    double dx;
    long double qx, w;

    /* initial value */
    dx = sin(pi * (2 * i - (n & 1) + 1.0) / (2 * n + 1.0));

    /* Newton iteration in double precision */
    for (k=0; ; k++) {
      double p0, p1, v;

      /* compute Pn(dx) */
      p0 = 1.0;
      p1 = dx;

      for (j=1; j < n; j++) {
	/* recurrence */
	v = ((2 * j + 1) * dx * p1 - j * p0) / (j + 1.0);

	p0 = p1;
	p1 = v;
      }

      /* v = Pn(dx) / (d Pn(dx) / d x) */
      v = p1 / (n * (p0 - dx * p1)) * (1.0 - dx * dx);

      /* check convergence */
      if (fabs(v) <= fabs(dx) * 1e-15)
	break;

      /* check overrun of iteration */
      if (k > 5) {
	fxt_error_set(FXT_ERROR_WARN,
		      "fxt_gauss_ll: Newton not converge %dd\n",
		      (int)i);
	break;
      }

      /* Newton update */
      dx -= v;
    }

    /* convert to long double */
    qx = (long double) dx;

    /* more Newton iteration in long double */
    for (k=0; ; k++) {
      long double p0, p1, v;

      /* compute Pn(x) */
      p0 = 1.0L;
      p1 = qx;

      for (j=1; j< n; j++) {
	/* recurrence */
	v = ((2 * j + 1) * qx * p1 - j * p0) / (j + 1.0L);

	p0 = p1;
	p1 = v;
      }

      /* v = Pn(x) / (d Pn(x) / d x) */
      v = p1 / (n * (p0 - qx * p1)) * (1.0L - qx * qx);

      /* check convergence */
      if (fabs((double)v) <= fabs((double)qx) * eps) {
	/* Gaussian weight */
	w = (1.0L - qx * qx) / (n * n * p0 * p0);

	break;
      }

      /* check overrun of iterations */
      if (k > 5) {
	fxt_error_set(FXT_ERROR_WARN,
		      "fxt_gauss_ll: Newton not converge %dq\n",
		      (int)i);

	/* Gaussian weight */
	w = (1.0L - qx * qx) / (n * n * p0 * p0);

	break;
      }

      /* Newton updat */
      qx -= v;
    }

    /* store Gauss nodes */
    xv->v[i + n - nhalf] = qx;
    xv->v[nhalf - 1 - i] = -qx;

    /* store Gauss weight */
    wv->v[i + n - nhalf] = w;
    wv->v[nhalf - 1 - i] = w;
  }

  /* set the results */
  *xp = xv;
  *wp = wv;
}
