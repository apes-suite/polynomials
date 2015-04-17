#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  pivoted LQ decomposition by Givens rotation
*********************************************************************/

#define MIN(x, y) ((x) < (y) ? (x) : (y))


/*** make givens rotation ***/
static double givens(double a, double b, double *c, double *s) {
  double z;

  if (b == 0.0) {
    *c = 1.0;
    *s = 0.0;
    z = a;
  } else if (fabs(b) > fabs(a)) {
    double t = a / b;
    double tt = sqrt(1.0 + t * t);
    *s = 1.0 / tt;
    *c = t * (*s);
    z = tt * b;
  } else {
    double t = b / a;
    double tt = sqrt(1.0 + t * t);
    *c = 1.0 / tt;
    *s = t * (*c);
    z = tt * a;
  }

  return z;
}


/*** encode givens rotation into a single value ***/
static double encode_givens(double c, double s) {
  double r;

  if (c == 0.0)
    r = 1.0;
  else if (fabs(s) < fabs(c))
    r = (c < 0 ? -s : s);
  else
    r = (s < 0 ? -1.0 : 1.0) / c;

  return r;
}


/*** decode givens rotation from a single value ***/
static void decode_givens(double r, double *c, double *s) {
  if (r == 1.0) {
    *c = 0.0;
    *s = 1.0;
  } else if (fabs(r) < 1.0) {
    *s = r;
    *c = sqrt(1.0 - r * r);
  } else {
    r = 1.0 / r;
    *c = r;
    *s = sqrt(1.0 - r * r);
  }
}


/*** LQ decomposition by Givens rotation ***/
void fxt_matld_glq(fxt_matld *mat, fxt_vecl *piv) {
  double *s;			/* norms of rows */
  double *ss;			/* arrays of sin */
  double *cc;			/* arrays of cos */
  long i, j, k;

  /* check null pointers */
  if (mat == NULL || piv == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq: null pointer\n");
    return;
  }

  /* check size */
  if (piv->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq: size mismatch\n");
    return;
  }

  /* allocate temporal arrays */
  s = (double *) malloc(sizeof(double) * mat->nrow);
  if (s == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_glq: allocation failed\n");
    return;
  }

  cc = (double *) malloc(sizeof(double) * mat->ncol);
  if (cc == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_glq: allocation failed\n");
    return;
  }

  ss = (double *) malloc(sizeof(double) * mat->ncol);
  if (ss == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_glq: allocation failed\n");
    return;
  }

  /* initializations */
  for (i=0; i< mat->nrow; i++) {
    piv->v[i] = i;

    /* squared 2-norm of i-th row */
    s[i] = 0.0;
    for (j=0; j< mat->ncol; j++)
      s[i] += MENT(mat, i, j) * MENT(mat, i, j);
  }

  /* Givens LQ decomposition */
  for (k=0; k< MIN(mat->nrow, mat->ncol); k++) {
    double t;
    long imax;

    /* pivot search */
    imax = k;
    for (i= k+1; i< mat->nrow; i++)
      if (s[imax] < s[i])
	imax = i;

    /* swap rows */
    if (imax != k) {
      t = s[imax];  s[imax] = s[k];  s[k] = t;

      fxt_matld_swaprow(mat, k, imax);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_vecl_swap(piv, k, imax);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    /* generate Givens rotation */
    t = MENT(mat, k, k);
    for (i = k+1; i < mat->ncol; i++) {
      /* compute Givens */
      t = givens(t, MENT(mat, k, i), &cc[i], &ss[i]);

      /* store encoded givens */
      MENT(mat, k, i) = encode_givens(cc[i], ss[i]);
    }
    MENT(mat, k, k) = t;

    /* apply Givens rotation */
    for (j = k+1; j < mat->nrow; j++) {

      t = MENT(mat, j, k);
      for (i= k+1; i < mat->ncol; i++) {
	double v = MENT(mat, j, i);
	MENT(mat, j, i) = v * cc[i] - t * ss[i];
	t = t * cc[i] + v * ss[i];
      }
      MENT(mat, j, k) = t;

      s[j] -= t * t;
    }
  }

  free(ss);
  free(cc);
  free(s);
}


/*** multiply a matrix by L^{-1} ***/
void fxt_matld_glq_mlinv(fxt_matld *mat, fxt_matld *lq) {
  long i, j, k;

  /* check null pointers */
  if (mat == NULL || lq == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq_mlinv: null pointer\n");
    return;
  }

  /* check size */
  if (mat->ncol != MIN(lq->nrow, lq->ncol)) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq_mlinv: size mismatch\n");
    return;
  }

  /* multiplication */
  for (k = mat->ncol - 1; k >= 0; k--) {
    /* diagonal of L */
    double p = 1.0 / MENT(lq, k, k);

    for (i=0; i< mat->nrow; i++) {
      /* 'solve' */
      double t = MENT(mat, i, k) *= p;

      for (j=0; j< k; j++)
	MENT(mat, i, j) -= t * MENT(lq, k, j);
    }
  }
}


/*** multiply a matrix by Q^{-1} ***/
void fxt_matld_glq_mqinv(fxt_matld *mat, fxt_matld *lq) {
  double *cc;			/* array of cos */
  double *ss;			/* array of sin */
  long i, j, k;

  /* check null pointers */
  if (mat == NULL || lq == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq_mqinv: null pointer\n");
    return;
  }

  /* check size */
  if (lq->ncol != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_glq_mqinv: size mismatch\n");
    return;
  }

  /* allocate memory */
  cc = (double *) malloc(sizeof(double) * lq->ncol);
  if (cc == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_glq_mqinv: allocation failed\n");
    return;
  }

  ss = (double *) malloc(sizeof(double) * lq->ncol);
  if (ss == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_glq_mqinv: allocation failed\n");
    return;
  }

  /* multiply Q^{-1} via Givens rotation */
  for (k=0; k< MIN(lq->ncol, lq->nrow); k++) {
    /* decode Givens rotation */
    for (i= k+1; i< lq->ncol; i++)
      decode_givens(MENT(lq, k, i), &cc[i], &ss[i]);

    /* apply Givens rotation */
    for (j=0; j< mat->nrow; j++) {
      double t = MENT(mat, j, k);

      for (i= k+1; i< lq->ncol; i++) {
	double v = MENT(mat, j, i);
	MENT(mat, j, i) = v * cc[i] - t * ss[i];
	t = t * cc[i] + v * ss[i];
      }

      MENT(mat, j, k) = t;
    }
  }

  /* deallocate temporal */
  free(ss);
  free(cc);
}
