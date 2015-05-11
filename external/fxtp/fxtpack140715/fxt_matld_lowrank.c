#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matld.h"

#include "fxt_vecl.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  low rank approximation of a matrix 
*********************************************************************/

static int lowrank_mgs(fxt_matld *a, double eps,
		       fxt_matld *x, fxt_matld *y);


/*** X * Y approximates A, returns rank ***/
int fxt_matld_lowrank(fxt_matld *a, double eps,
		      fxt_matld *x, fxt_matld *y) {

  /* check null pointers */
  if (a == NULL || x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lowrank: null pointer\n");
    return 0;
  }

  return lowrank_mgs(a, eps, x, y);
}


static double threshould(fxt_matld *a, double eps);

/*** low-rank approximation by Modified Gramm-Schmidt ***/
static int lowrank_mgs(fxt_matld *a, double eps,
			fxt_matld *x, fxt_matld *r) {
  fxt_matld *q;			/* the approximation */
  fxt_vecl *idx;		/* ordering */
  long i, j, k, kmax;
  double s, ssum;
  int retc = 0;

  /* check size */
  if (x->ncol != r->nrow ||
      a->nrow != x->nrow || a->ncol != r->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lowrank: size mismatch\n");
    return retc;
  }

  /* the maximum rank assumed */
  kmax = x->ncol;

  /* initialize QR decomposition */
  q = fxt_matld_clone(a, 0, a->nrow - 1, 0, a->ncol - 1);
  if (q == NULL) {
    fxt_error_raise();
    return retc;
  }

  /* the permutation array */
  idx = fxt_vecl_new(a->ncol);
  if (idx == NULL) {
    fxt_error_raise();
    return retc;
  }

  /* initilize permutation */
  for (i=0; i< idx->n; i++)
    idx->v[i] = i;

  /* the Gramm-Schmidt loop */
  for (k = 0; ; k ++) {
    long imax = k;		/* pivot index */
    double smax = 0.0;		/* norm of pivot coloumn */
    ssum = 0.0;			/* total Frobenius norm */

    /* pivot search */
    for (i = k; i < a->ncol; i++) {
      /* norm of this column */
      s = 0.0;
      for (j = 0; j < a->nrow; j++)
	s += MENT(q, j, i) * MENT(q, j, i);

      /* check maximality */
      if (i == k || smax < s) {
	smax = s;
	imax = i;
      }

      /* add up norm */
      ssum += s;
    }

    /* check for convergence */
    if (ssum <= eps * eps)
      break;

    /* check overrun of the loop */
    if (k >= kmax) {
      retc = 1;			/* insufficient size k */
      break;
    }

    /* swap pivot column */
    if (imax != k) {
      /* swap matrix column */
      fxt_matld_swapcol(q, k, imax);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return retc;

      /* update index */
      fxt_vecl_swap(idx, k, imax);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return retc;
    }

    /* scale pivot column */
    s = 1.0 / sqrt(smax);
    for (j = 0; j < a->nrow; j++)
      MENT(q, j, k) *= s;

    /* clear lower triangular part of R */
    for (i = 0; i < k; i++)
      MENT(r, k, idx->v[i]) = 0.0;

    /* diagonal part of R */
    MENT(r, k, idx->v[k]) = sqrt(smax);

    /* orthogonalization */
    for (i = k + 1; i < a->ncol; i++) {
      /* inner product */
      s = 0.0;
      for (j=0; j< a->nrow; j++)
	s += MENT(q, j, k) * MENT(q, j, i);

      /* upper triangular part of R */
      MENT(r, k, idx->v[i]) = s;

      /* orthogonalize columns */
      for (j=0; j< a->nrow; j++)
	MENT(q, j, i) -= s * MENT(q, j, k);
    }
  }

  /* resize matrices */
  fxt_matld_shrink(x, a->nrow, k);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return retc;

  fxt_matld_shrink(r, k, a->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return retc;

  /* copy Q into x */
  if (k > 0) {
    fxt_matld_setprt(x, q, 0, a->nrow - 1, 0, k - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return retc;
  }

  /* drop small entries */
  s = eps * eps - ssum;
  if (s > 0.0) {
    double th = threshould(r, sqrt(s));
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return retc;

    fxt_matld_drop(r, th);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return retc;
  }

  /* deallocate temporals */
  fxt_vecl_del(idx);
  fxt_matld_del(q);

  /* return the rank */
  return retc;
}


/*** drop threshould in abs-sum ***/
double threshould(fxt_matld *mat, double eps) {
  double th;			/* the threshould */
  int c;			/* loop counter */

  /* initial guess */
  th = eps;

  for (c=0; ; c++) {
    long i;
    double toterr = 0.0;	/* total error by dropping */
    double maxerr = 0.0;	/* maximum error by dropping */

    /* try dropping */
    for (i = 0; i < mat->nrow * mat->ncol; i++) {
      double v = fabs(mat->v[i]);

      if (v <= th) {		/* dropped entry */
	/* update total error */
	toterr += v;

	/* update maximum error */
	if (maxerr < v)
	  maxerr = v;
      }
    }

    /* check convergence */
    if (toterr < eps)
      break;

    /* update threshould, linear extrapolation */
    if (eps / toterr > 0.7)
      th = th * 0.7;
    else
      th = th * eps / toterr;

    /* modify if it is too large */
    if (th > maxerr)
      th = maxerr * (1 - 1e-15);
  }

  return th;
}
