#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fxt_error.h"
#include "fxt_vecl.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for array of long int
*********************************************************************/

/*** create a new array ***/
fxt_vecl* fxt_vecl_new(long size) {
  fxt_vecl* x;			/* the vector */

  /* check input */
  if (size < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_new: negative size\n");
    return NULL;
  }

  /* allocate memory for the structure */
  x = (fxt_vecl*) malloc(sizeof(fxt_vecl));
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecl_new: allocation failed\n");
    return NULL;
  }

  /* allocate memory for the vector */
  x->n = size;
  if (size == 0)
    x->v = NULL;
  else {
    x->v = (long*) malloc(sizeof(long) * size);
    if (x->v == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_vecl_new: allocation failed\n");
      return NULL;
    }
  }

  return x;
}


/*** deallocate memory ***/
void fxt_vecl_del(fxt_vecl *x) {

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_del: null pointer\n");
    return;
  }

  /* deallocate vector */
  if (x->v != NULL) {
    free(x->v);
    x->v = NULL;
  }

  /* deallocate structure */
  free(x);
}


/*** create a clone of specified range ***/
fxt_vecl* fxt_vecl_clone(fxt_vecl *x, long ini, long fin) {
  fxt_vecl *y;			/* the clone */
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_clone: null pointer\n");
    return NULL;
  }

  /* check input */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_clone: irregal range\n");
    return NULL;
  }

  /* create the new vector */
  y = fxt_vecl_new(fin - ini + 1);
  if (y == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i= ini; i <= fin; i++)
    y->v[i - ini] = x->v[i];

  return y;
}


/*** copy y into a range of x ***/
void fxt_vecl_prtset(fxt_vecl *x, long ini, long fin, fxt_vecl *y) {
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_prtset: null pointer\n");
    return;
  }

  /* check input */
  if (ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_set: irregal range\n");
    return;
  }

  if (y->n != fin - ini + 1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_set: size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< y->n; i++)
    x->v[i + ini] = y->v[i];
}


/*** swap two entries ***/
void fxt_vecl_swap(fxt_vecl *x, long i, long j) {
  long t;			/* temporal */

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_swap: null pointer\n");
    return;
  }

  /* check input */
  if (i < 0 || x->n <= i || j < 0 || x->n <= j) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_swap: irregal index(s)\n");
    return;
  }

  /* swap */
  t = x->v[i];
  x->v[i] = x->v[j];
  x->v[j] = t;
}


/*** permute the vector ***/
void fxt_vecl_perm(fxt_vecl *x, fxt_vecl *perm) {
  long *idx;			/* original index of the entry */
  long i, ii;

  /* check null pointers */
  if (x == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_perm: null pointer\n");
    return;
  }

  /* allocate index array */
  idx = (long *) malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecl_perm: allocation failed\n");
    return;
  }

  /* check input */
  if (perm->n != x->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecl_perm: size mismatch\n");
    return;
  }

  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i] ||
	(++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecl_perm: error in permutation vector\n");
      return;
    }

  /* initialize index array */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute vector in the the same way as for a matrix */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the value didn't come yet */
      /* save the current value */
      long t = x->v[ii];

      i = ii;			/* index to fix */
      while (perm->v[i] != ii) {
	long k = perm->v[i];

	/* get the value */
	x->v[i] = x->v[k];

	/* update the index */
	idx[i] = k;

	/* goto next position */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* get the value */
      x->v[i] = t;

      /* update the index */
      idx[i] = ii;
    }

  /* deallocate temporal */
  free(idx);
}
