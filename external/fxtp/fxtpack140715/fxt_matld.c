#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_matld.h"

#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for matrix of double, indexed with long int
*********************************************************************/

#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*** create a new matrix ***/
fxt_matld* fxt_matld_new(long nrow, long ncol) {
  fxt_matld *mat;

  /* check size */
  if (nrow < 0 || ncol < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_new: negative size\n");
    return NULL;
  }

  /* allocate data structure */
  mat = (fxt_matld*) malloc(sizeof(fxt_matld));
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_new: allocation failed\n");
    return NULL;
  }

  /* set sizes */
  mat->nrow = nrow;
  mat->ncol = ncol;

  /* allocate array */
  if (nrow == 0 || ncol == 0)
    mat->v = NULL;
  else {
    mat->v = (double*) malloc(sizeof(double) * nrow * ncol);
    if (mat->v == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_matld_new: allocation failed\n");
      return NULL;
    }
  }

  /* return data */
  return mat;
}


/*** shrink matrix ***/
void fxt_matld_shrink(fxt_matld *mat, long nrow, long ncol) {
  
  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_shrink: null pointer\n");
    return;
  }

  /* check size */
  if (nrow > mat->nrow || ncol > mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_shrink: larger size required\n");
    return;
  }

  if (0 < nrow && 0 < ncol &&
      (nrow < mat->nrow || ncol < mat->ncol)) {
    long i, j;

    /* copy old value */
    for (i=0; i< nrow; i++)
      for (j=0; j< ncol; j++)
	mat->v[i * ncol + j] = MENT(mat, i, j);
  }

  /* resize */
  mat->nrow = nrow;
  mat->ncol = ncol;
}


/*** deallocate memory ***/
void fxt_matld_del(fxt_matld *mat) {

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_del: null pointer\n");
    return;
  }

  /* deallocate array */
  if (mat->v != NULL) {
    free(mat->v);
    mat->v = NULL;
  }

  /* deallocate data structure */
  free(mat);
}


/*** create a clone of the range ***/
fxt_matld* fxt_matld_clone(fxt_matld *mato, long rini, long rfin,
			   long cini, long cfin) {
  fxt_matld *matc;		/* the clone */
  long i, j;

  /* check null pointers */
  if (mato == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_clone: null pointer\n");
    return NULL;
  }

  /* check range */
  if (rfin < rini || rini < 0 || mato->nrow <= rfin ||
      cfin < cini || cini < 0 || mato->ncol <= cfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_clone: range error\n");
    return NULL;
  }

  /* create a new matrix */
  matc = fxt_matld_new(rfin - rini + 1, cfin - cini + 1);
  if (matc == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i= rini; i <= rfin; i++)
    for (j= cini; j <= cfin; j++)
      MENT(matc, i - rini, j - cini) = MENT(mato, i, j);

  /* return result */
  return matc;
}


/*** fill zero to a submatrix */
void fxt_matld_zero(fxt_matld *mat, long rini, long rfin,
		    long cini, long cfin) {
  long i, j;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_zero: null pointer\n");
    return;
  }

  /* check range */
  if (rfin < rini || rini < 0 || mat->nrow <= rfin ||
      cfin < cini || cini < 0 || mat->ncol <= cfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_zero: range error\n");
    return;
  }

  /* fill zeros */
  for (i= rini; i <= rfin; i++)
    for (j= cini; j <= cfin; j++)
      MENT(mat, i, j) = 0.0;
}


/*** matrix copy: A[r0:r1, c0:c1] = B ***/
void fxt_matld_prtset(fxt_matld *a,
		      long r0, long r1, long c0, long c1,
		      fxt_matld *b) {
  long i, j;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtset: null pointer\n");
    return;
  }

  /* check range */
  if (r1 < r0 || r0 < 0 || a->nrow <= r1 ||
      c1 < c0 || c0 < 0 || a->ncol <= c1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtset: irregal range\n");
    return;
  }

  /* check size */
  if (r1 - r0 + 1 != b->nrow || c1 - c0 + 1 != b->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtset: size mismatch\n");
    return;
  }

  /* copy values */
  for (i= r0; i<= r1; i++)
    for (j= c0; j<= c1; j++)
      MENT(a, i, j) = MENT(b, i-r0, j-c0);
}


/*** matrix copy: A = B[r0:r1, c0:c1] ***/
void fxt_matld_setprt(fxt_matld *a, fxt_matld *b,
		      long r0, long r1, long c0, long c1) {
  long i, j;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setprt: null pointer\n");
    return;
  }

  /* check range */
  if (r1 < r0 || r0 < 0 || b->nrow <= r1 ||
      c1 < c0 || c0 < 0 || b->ncol <= c1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setprt: irregal range\n");
    return;
  }

  /* check size */
  if (a->nrow != r1 - r0 + 1 || a->ncol != c1 - c0 + 1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setprt: size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< a->nrow; i++)
    for (j=0; j< a->ncol; j++)
      MENT(a, i, j) = MENT(b, i + r0, j + c0);
}


/*** matrix copy: A[r0:r1, c0:c1] = B[r2:r3, c2:c3] ***/
void fxt_matld_prtsetprt(fxt_matld *a,
			 long r0, long r1, long c0, long c1,
			 fxt_matld *b,
			 long r2, long r3, long c2, long c3){
  long i, j;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetprt: null pointer\n");
    return;
  }

  /* check range */
  if (r1 < r0 || r0 < 0 || a->nrow <= r1 ||
      c1 < c0 || c0 < 0 || a->ncol <= c1 ||
      r3 < r2 || r2 < 0 || b->nrow <= r3 ||
      c3 < c2 || c2 < 0 || b->ncol <= c3) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetprt: irregal range\n");
    return;
  }

  /* check size */
  if (r1 - r0 != r3 - r2 || c1 - c0 != c3 - c2) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetprt: size mismatch\n");
    return;
  }

  /* compute difference */
  r2 -= r0;  c2 -= c0;

  /* copy values */
  for (i= r0; i <= r1; i++)
    for (j= c0; j <= c1; j++)
      MENT(a, i, j) = MENT(b, i+r2, j+c2);
}


/*** set to identity matrix ***/
void fxt_matld_unit(fxt_matld *mat) {
  long i, j;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_unit: null pointer\n");
    return;
  }

  /* check size */
  if (mat->nrow != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_unit: non-square matrix\n");
    return;
  }

  /* fill values */
  for (i=0; i< mat->nrow; i++) {
    /* first fill zero */
    for (j=0; j< mat->ncol; j++)
      MENT(mat, i, j) = 0.0;

    /* make diagonal one */
    MENT(mat, i, i) = 1.0;
  }
}


/*** matrix subtraction a - b ***/
void fxt_matld_sub(fxt_matld *a, fxt_matld *b) {
  long i;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_sub: null pointer\n");
    return;
  }

  /* check size */
  if (a->nrow != b->nrow || a->ncol != b->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_sub: size mismatch\n");
    return;
  }

  /* subtraction */
  for (i=0; i< a->nrow * a->ncol; i++)
    a->v[i] -= b->v[i];
}


/*** matrix-matrix multiplication c = a * b ***/
void fxt_matld_mul(fxt_matld *c, fxt_matld *a, fxt_matld *b) {
  long nrow, ncol, nmid;

  /* check null pointers */
  if (a == NULL || b == NULL || c == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_mul: null pointer\n");
    return;
  }

  /* check size */
  if (c->nrow != a->nrow || a->ncol != b->nrow ||
      b->ncol != c->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_mul: size mismatch\n");
    return;
  }

  nrow = c->nrow;
  ncol = c->ncol;
  nmid = a->ncol;

  /* multiplication */
  if (nmid == 0) {		/* zero-sized product */
    if (nrow > 0 && ncol > 0) {
      fxt_matld_zero(c, 0, nrow - 1, 0, ncol - 1);
      fxt_error_raise();
    }

  } else {			/* usual case */
    long i, j, k;

    for (i=0; i < nrow; i++) {
      /* the first term */
      for (j=0; j < ncol; j++)
	MENT(c, i, j) = MENT(a, i, 0) * MENT(b, 0, j);

      /* the other terms */
      for (k=1; k < nmid; k++)
	for (j=0; j< ncol; j++)
	  MENT(c, i, j) += MENT(a, i, k) * MENT(b, k, j);
    }
  }
}


/*** Frobenius norm ***/
double fxt_matld_normf(fxt_matld *mat) {
  double s;
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_normf: null pointer\n");
    return 0.0;
  }

  /* check zero-sized matrix */
  if (mat->nrow == 0 || mat->ncol == 0)
    return 0.0;

  /* compute sum of square */
  s = mat->v[0] * mat->v[0];
  for (i=1; i< mat->nrow * mat->ncol; i++)
    s += mat->v[i] * mat->v[i];

  /* return norm */
  return sqrt(s);
}


/*** 2-norm by power method ***/
double fxt_matld_norm2(fxt_matld *mat, double prec) {
  fxt_vecld *u, *v;		/* eigen vectors */
  double s, s0, ss;		/* squared norm */
  int c, cc;			/* iteration counter */

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2: null pointer\n");
    return 0.0;
  }

  /* check size */
  if (mat->nrow == 0 || mat->ncol == 0)
    return 0.0;

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2: error in request precision\n");
    return 0.0;
  }

  /* allocate eigenvectors */
  u = fxt_vecld_new(mat->nrow);
  if (u == NULL) {
    fxt_error_raise();
    return 0.0;
  }

  v = fxt_vecld_new(mat->ncol);
  if (v == NULL) {
    fxt_error_raise();
    return 0.0;
  }

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    do {

      /* initialize by random vector */
      fxt_vecld_rand(v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

    } while (fxt_vecld_norm2(v) == 0.0);

    /* initialize norm */
    s = s0 = 0.0;

    /* power method iteration */
    for (c = 0; ; c ++) {
      fxt_vecld_normalize(v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* v = A^T A v */
      fxt_matld_setmatvec(u, mat, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_settrmatvec(v, mat, u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* get the norm */
      s = fxt_vecld_norm2(v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      /* check convergence */
      if (s == 0.0 || (c > MINPOWERITER && (s - s0) / s <= prec))
	break;

      s0 = s;
    }

    s = MAX(s, s0);

    if (s == 0.0 || fabs(s - ss) / MAX(s, ss) <= prec)
      break;

    ss = MAX(s, ss);
  }

  s = MAX(s, ss);

  /* deallocate temporals */
  fxt_vecld_del(v);
  fxt_vecld_del(u);

  /* return norm */
  return sqrt(s);
}


/*** swap rows ***/
void fxt_matld_swaprow(fxt_matld *mat, long j, long k) {
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_swaprow: null pointer\n");
    return;
  }

  /* check range */
  if (j < 0 || mat->nrow <= j ||
      k < 0 || mat->nrow <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_swaprow: index out of range\n");
    return;
  }

  /* check no swap */
  if (j == k)
    return;

  /* swap rows */
  for (i=0; i< mat->ncol; i++) {
    double t = MENT(mat, j, i);

    MENT(mat, j, i) = MENT(mat, k, i);

    MENT(mat, k, i) = t;
  }
}


/*** swap columns ***/
void fxt_matld_swapcol(fxt_matld *mat, long j, long k) {
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_swapcol: null pointer\n");
    return;
  }

  /* check range */
  if (j < 0 || mat->ncol <= j ||
      k < 0 || mat->ncol <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_swapcol: index out of range\n");
    return;
  }

  /* check no swap */
  if (j == k)
    return;

  /* swap columns */
  for (i=0; i< mat->nrow; i++) {
    double t = MENT(mat, i, j);

    MENT(mat, i, j) = MENT(mat, i, k);

    MENT(mat, i, k) = t;
  }
}


/*** permute rows ***/
void fxt_matld_permrow(fxt_matld *mat, fxt_vecl *perm) {
  long i, j, ii;
  long *idx;			/* original index */
  double *t;			/* temporal row */

  /* check null pointers */
  if (mat == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_permrow: null pointer\n");
    return;
  }

  /* check size */
  if (perm->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_permrow: size mismatch\n");
    return;
  }

  /* check range */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_permrow: error in perm vector\n");
      return;
    }

  /* check no permutation */
  if (perm->n <= 1 || mat->ncol == 0)
    return;

  /* allocate index */
  idx = (long*) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_permrow: allocation failed\n");
    return;
  }

  /* allocate temporal row */
  t = (double*) malloc(sizeof(double) * mat->ncol);
  if (t == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_permrow: allocation failed\n");
    return;
  }

  /* check permutation */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_permrow: error in perm vector\n");
      return;
    }

  /* initialize index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the row didn't come yet */
      /* save the current row */
      for (j=0; j< mat->ncol; j++)
	t[j] = MENT(mat, ii, j);

      i = ii;
      while (perm->v[i] != ii) {
	/* the target index */
	long k = perm->v[i];

	/* get the row */
	for (j=0; j< mat->ncol; j++)
	  MENT(mat, i, j) = MENT(mat, k, j);

	/* update index */
	idx[i] = k;

	/* next index */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* get the row */
      for (j=0; j< mat->ncol; j++)
	MENT(mat, i, j) = t[j];

      /* update index */
      idx[i] = ii;
    }

  /* deallocate temporal arrays */
  free(t);
  free(idx);
}


/*** permute columns ***/
void fxt_matld_permcol(fxt_matld *mat, fxt_vecl *perm) {
  long i, j, ii;
  long *idx;			/* original index */
  double *t;			/* temporal column */

  /* check null pointers */
  if (mat == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_permcol: null pointer\n");
    return;
  }

  /* check size */
  if (perm->n != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_permcol: size mismatch\n");
    return;
  }

  /* check range */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_permcol: error in perm vector\n");
      return;
    }

  /* check no permutation */
  if (perm->n <= 1 || mat->nrow == 0)
    return;

  /* allocate index array */
  idx = (long*) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_permcol: allocation failed\n");
    return;
  }

  /* allocate temporal column */
  t = (double*) malloc(sizeof(double) * mat->nrow);
  if (t == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_permcol: allocation failed\n");
    return;
  }

  /* check permutation */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_permcol: error in perm vector\n");
      return;
    }

  /* initialize index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the column didn't come yet */
      /* save the current column */
      for (j=0; j< mat->nrow; j++)
	t[j] = MENT(mat, j, ii);

      i = ii;
      while (perm->v[i] != ii) {
	/* the target index */
	long k = perm->v[i];

	/* get the column */
	for (j=0; j< mat->nrow; j++)
	  MENT(mat, j, i) = MENT(mat, j, k);

	/* update index */
	idx[i] = k;

	/* next index */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* get the column */
      for (j=0; j< mat->nrow; j++)
	MENT(mat, j, i) = t[j];

      /* update index */
      idx[i] = ii;
    }

  /* deallocate temporal arrays */
  free(t);
  free(idx);
}


/*** inverse permutation of the rows ***/
void fxt_matld_ipermrow(fxt_matld *mat, fxt_vecl *perm) {
  long i, j, ii;
  long *idx;			/* original index */
  double *t;			/* temporal row */

  /* check null pointers */
  if (mat == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_ipermrow: null pointer\n");
    return;
  }

  /* check size */
  if (perm->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_ipermrow: size mismatch\n");
    return;
  }

  /* check range */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_ipermrow: error in perm vector\n");
      return;
    }

  /* check no permutation */
  if (perm->n <= 1 || mat->ncol == 0)
    return;

  /* allocate index */
  idx = (long*) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_ipermrow: allocation failed\n");
    return;
  }

  /* allocate temporal row */
  t = (double*) malloc(sizeof(double) * mat->ncol);
  if (t == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_ipermrow: allocation failed\n");
    return;
  }

  /* check permutation */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_ipermrow: error in perm vector\n");
      return;
    }

  /* initialize index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the row didn't moved yet */
      /* save the current row */
      for (j=0; j< mat->ncol; j++)
	t[j] = MENT(mat, ii, j);

      i = ii;
      while (perm->v[i] != ii) {
	/* the target index */
	long k = perm->v[i];

	/* put the row */
	for (j=0; j< mat->ncol; j++) {
	  double tt = t[j];
	  t[j] = MENT(mat, k, j);
	  MENT(mat, k, j) = tt;
	}

	/* update index */
	idx[i] = k;

	/* next index */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* put the row */
      for (j=0; j< mat->ncol; j++)
	MENT(mat, ii, j) = t[j];

      /* update index */
      idx[i] = ii;
    }

  /* deallocate temporal arrays */
  free(t);
  free(idx);
}


/*** get k-th row ***/
void fxt_matld_getrow(fxt_vecld *v, fxt_matld *mat, long k) {
  long i;

  /* check null pointers */
  if (v == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getrow: null pointer\n");
    return;
  }

  /* check input */
  if (k < 0 || mat->nrow <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getrow: index out of range\n");
    return;
  }

  /* check size */
  if (v->n != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getrow: vector size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< mat->ncol; i++)
    v->v[i] = MENT(mat, k, i);
}


/*** get k-th column ***/
void fxt_matld_getcol(fxt_vecld *v, fxt_matld *mat, long k) {
  long i;

  /* check null pointers */
  if (v == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getcol: null pointer\n");
    return;
  }

  /* check input */
  if (k < 0 || mat->ncol <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getcol: index out of range\n");
    return;
  }

  /* check size */
  if (v->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_getcol: vector size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< mat->nrow; i++)
    v->v[i] = MENT(mat, i, k);
}


/*** scale rows (multiply diagonal matrix from right) ***/
void fxt_matld_scalerow(fxt_matld *mat, fxt_vecld *sc) {
  long i, j;

  /* check null pointers */
  if (mat == NULL || sc == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_scalerow: null pointer\n");
    return;
  }

  /* check size */
  if (sc->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_scalerow: vector size mismatch\n");
    return;
  }

  /* scale values */
  for (i=0; i< mat->nrow; i++)
    for (j=0; j< mat->ncol; j++)
      MENT(mat, i, j) *= sc->v[i];
}


/*** scale columns (multiply diagonal matrix from left) ***/
void fxt_matld_scalecol(fxt_matld *mat, fxt_vecld *sc) {
  long i, j;

  /* check null pointers */
  if (mat == NULL || sc == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_scalecol: null pointer\n");
    return;
  }

  /* check size */
  if (sc->n != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_scalecol: vector size mismatch\n");
    return;
  }

  /* scale values */
  for (i=0; i< mat->nrow; i++)
    for (j=0; j< mat->ncol; j++)
      MENT(mat, i, j) *= sc->v[j];
}


/*** 2-norm of the rows ***/
void fxt_matld_norm2rows(fxt_vecld *rnorm, fxt_matld *mat) {
  long i, j;

  /* check null pointers */
  if (rnorm == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2rows: null pointer\n");
    return;
  }

  /* check size */
  if (rnorm->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2rows: size mismatch\n");
    return;
  }

  /* get norms */
  for (i=0; i< mat->nrow; i++) {
    double s = 0.0;

    for (j=0; j< mat->ncol; j++)
      s += MENT(mat, i, j) * MENT(mat, i, j);

    rnorm->v[i] = sqrt(s);
  }
}


/*** 2-norm of the columns ***/
void fxt_matld_norm2cols(fxt_vecld *cnorm, fxt_matld *mat) {
  long i, j;

  /* check null pointers */
  if (cnorm == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2cols: null pointer\n");
    return;
  }

  /* check size */
  if (cnorm->n != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_norm2cols: size mismatch\n");
    return;
  }

  /* zeroing vector */
  fxt_vecld_zero(cnorm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* get norms */
  for (i=0; i< mat->nrow; i++)
    for (j=0; j< mat->ncol; j++)
      cnorm->v[j] += MENT(mat, i, j) * MENT(mat, i, j);

  /* square root */
  fxt_vecld_sqrt(cnorm);
  fxt_error_raise();
}
