#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matll.h"

#include "fxt_vecll.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for matrix of long double, indexed with long int
*********************************************************************/

/*** create a new matrix ***/
fxt_matll* fxt_matll_new(long nrow, long ncol) {
  fxt_matll *mat;

  /* check size */
  if (nrow < 0 || ncol < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_new: negative size\n");
    return NULL;
  }

  /* allocate data structure */
  mat = (fxt_matll*) malloc(sizeof(fxt_matll));
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matll_new: allocation failed\n");
    return NULL;
  }

  /* set sizes */
  mat->nrow = nrow;
  mat->ncol = ncol;

  /* allocate array */
  if (nrow == 0 || ncol == 0)
    mat->v = NULL;
  else {
    mat->v = (long double*) malloc(sizeof(long double) * nrow * ncol);
    if (mat->v == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_matll_new: allocation failed\n");
      return NULL;
    }
  }

  /* return data */
  return mat;
}


/*** deallocate memory ***/
void fxt_matll_del(fxt_matll *mat) {

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_del: null pointer\n");
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
fxt_matll* fxt_matll_clone(fxt_matll *mato, long rini, long rfin,
			   long cini, long cfin) {
  fxt_matll *matc;		/* the clone */
  long i, j;

  /* check null pointers */
  if (mato == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_del: null pointer\n");
    return NULL;
  }

  /* check range */
  if (rfin < rini || rini < 0 || mato->nrow <= rfin ||
      cfin < cini || cini < 0 || mato->ncol <= cfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_clone: range error\n");
    return NULL;
  }

  /* create a new matrix */
  matc = fxt_matll_new(rfin - rini + 1, cfin - cini + 1);
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


/*** matrix copy: A[r0:r1, c0:c1] = B ***/
void fxt_matll_set(fxt_matll *a, long r0, long r1, long c0, long c1,
		   fxt_matll *b) {
  long i, j;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_set: null pointer\n");
    return;
  }

  /* check range */
  if (r1 < r0 || r0 < 0 || a->nrow <= r1 ||
      c1 < c0 || c0 < 0 || a->ncol <= c1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_set: irregal range\n");
    return;
  }

  /* check size */
  if (r1 - r0 + 1 != b->nrow || c1 - c0 + 1 != b->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_set: size mismatch\n");
    return;
  }

  /* copy values */
  for (i= r0; i<= r1; i++)
    for (j= c0; j<= c1; j++)
      MENT(a, i, j) = MENT(b, i-r0, j-c0);
}


/*** matrix subtraction a - b ***/
void fxt_matll_sub(fxt_matll *a, fxt_matll *b) {
  long i;

  /* check null pointers */
  if (a == NULL || b == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_del: null pointer\n");
    return;
  }

  /* check size */
  if (a->nrow != b->nrow || a->ncol != b->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_sub: size mismatch\n");
    return;
  }

  /* subtraction */
  for (i=0; i< a->nrow * a->ncol; i++)
    a->v[i] -= b->v[i];
}


/*** matrix-matrix multiplication c = a * b ***/
void fxt_matll_mul(fxt_matll *c, fxt_matll *a, fxt_matll *b) {
  long nrow, ncol, nmid;

  /* check null pointers */
  if (a == NULL || b == NULL || c == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_mul: null pointer\n");
    return;
  }

  /* check size */
  if (c->nrow != a->nrow || a->ncol != b->nrow ||
      b->ncol != c->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_mul: size mismatch\n");
    return;
  }

  nrow = c->nrow;
  ncol = c->ncol;
  nmid = a->ncol;

  /* multiplication */
  if (nmid == 0) {		/* zero-sized product */
    long i, j;

    if (nrow > 0 && ncol > 0)
      for (i=0; i< nrow; i++)
	for (j=0; j< ncol; j++)
	  MENT(c, i, j) = 0.0;

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
long double fxt_matll_normf(fxt_matll *mat) {
  long double s;
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_normf: null pointer\n");
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


/*** permute rows ***/
void fxt_matll_permrow(fxt_matll *mat, fxt_vecl *perm) {
  long i, j, ii;
  long *idx;			/* original index */
  long double *t;			/* temporal row */

  /* check null pointers */
  if (mat == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_permrow: null pointer\n");
    return;
  }

  /* check size */
  if (perm->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_permrow: size mismatch\n");
    return;
  }

  /* check range */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matll_permrow: error in perm vector\n");
      return;
    }

  /* check no permutation */
  if (perm->n <= 1 || mat->ncol == 0)
    return;

  /* allocate index */
  idx = (long*) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matll_permrow: allocation failed\n");
    return;
  }

  /* allocate temporal row */
  t = (long double*) malloc(sizeof(long double) * mat->ncol);
  if (t == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matll_permrow: allocation failed\n");
    return;
  }

  /* check permutation */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matll_permrow: error in perm vector\n");
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


/*** inverse permutation of the rows ***/
void fxt_matll_ipermrow(fxt_matll *mat, fxt_vecl *perm) {
  long i, j, ii;
  long *idx;			/* original index */
  long double *t;			/* temporal row */

  /* check null pointers */
  if (mat == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_ipermrow: null pointer\n");
    return;
  }

  /* check size */
  if (perm->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_ipermrow: size mismatch\n");
    return;
  }

  /* check range */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matll_ipermrow: error in perm vector\n");
      return;
    }

  /* check no permutation */
  if (perm->n <= 1 || mat->ncol == 0)
    return;

  /* allocate index */
  idx = (long*) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matll_ipermrow: allocation failed\n");
    return;
  }

  /* allocate temporal row */
  t = (long double*) malloc(sizeof(long double) * mat->ncol);
  if (t == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matll_ipermrow: allocation failed\n");
    return;
  }

  /* check permutation */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matll_ipermrow: error in perm vector\n");
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
	  long double tt = t[j];
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
void fxt_matll_getrow(fxt_vecll *v, fxt_matll *mat, long k) {
  long i;

  /* check null pointers */
  if (v == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getrow: null pointer\n");
    return;
  }

  /* check input */
  if (k < 0 || mat->nrow <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getrow: index out of range\n");
    return;
  }

  /* check size */
  if (v->n != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getrow: vector size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< mat->ncol; i++)
    v->v[i] = MENT(mat, k, i);
}


/*** get k-th column ***/
void fxt_matll_getcol(fxt_vecll *v, fxt_matll *mat, long k) {
  long i;

  /* check null pointers */
  if (v == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getcol: null pointer\n");
    return;
  }

  /* check input */
  if (k < 0 || mat->ncol <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getcol: index out of range\n");
    return;
  }

  /* check size */
  if (v->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_getcol: vector size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< mat->nrow; i++)
    v->v[i] = MENT(mat, i, k);
}


/*** scale rows (multiply diagonal matrix from right) ***/
void fxt_matll_scalerow(fxt_matll *mat, fxt_vecll *sc) {
  long i, j;

  /* check null pointers */
  if (mat == NULL || sc == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_scalerow: null pointer\n");
    return;
  }

  /* check size */
  if (sc->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_scalerow: vector size mismatch\n");
    return;
  }

  /* scale values */
  for (i=0; i< mat->nrow; i++)
    for (j=0; j< mat->ncol; j++)
      MENT(mat, i, j) *= sc->v[i];
}


/*** conversion to double ***/
fxt_matld* fxt_matll_to_matld(fxt_matll *a) {
  fxt_matld *b;
  long i;

  /* check null pointers */
  if (a == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_to_matld: null pointer\n");
    return NULL;
  }

  /* generate new matrix */
  b = fxt_matld_new(a->nrow, a->ncol);
  if (b == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i=0; i< a->nrow * a->ncol; i++)
    b->v[i] = (double) a->v[i];

  /* return result */
  return b;
}


/*** conversion from double ***/
fxt_matll* fxt_matld_to_matll(fxt_matld *a) {
  fxt_matll *b;
  long i;

  /* check null pointers */
  if (a == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_to_matll: null pointer\n");
    return NULL;
  }

  /* generate new matrix */
  b = fxt_matll_new(a->nrow, a->ncol);
  if (b == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i=0; i< a->nrow * a->ncol; i++)
    b->v[i] = (long double) a->v[i];

  /* return result */
  return b;
}
