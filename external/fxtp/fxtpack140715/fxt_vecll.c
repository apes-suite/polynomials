#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_vecll.h"

#include "fxt_math.h"
#include "fxt_vecl.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for vector of long double, indexed with long int
*********************************************************************/

/*** create a new vector ***/
fxt_vecll* fxt_vecll_new(long size) {
  fxt_vecll *x;

  /* check input */
  if (size < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_new: negative size\n");
    return NULL;
  }

  /* allocate data structure */
  x = (fxt_vecll*) malloc(sizeof(fxt_vecll));
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecll_new: allocation failed\n");
    return NULL;
  }

  /* set size */
  x->n = size;

  /* allocate vector */
  if (size == 0)
    x->v = NULL;
  else {
    x->v = (long double*) malloc(sizeof(long double) * x->n);
    if (x->v == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_vecll_new: allocation failed\n");
      return NULL;
    }
  }

  /* return data */
  return x;
}


/*** deallocate memory ***/
void fxt_vecll_del(fxt_vecll *x) {

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_del: null pointer\n");
    return;
  }

  /* deallocate vector */
  if (x->v != NULL) {
    free(x->v);
    x->v = NULL;
  }

  /* deallocate data structure */
  free(x);
}


/*** get a clone ***/
fxt_vecll* fxt_vecll_clone(fxt_vecll *x, long ini, long fin) {
  fxt_vecll *y;
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_clone: null pointer\n");
    return NULL;
  }

  /* check input */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_clone: irregal range\n");
    return NULL;
  }

  /* get new vector */
  y = fxt_vecll_new(fin - ini + 1);
  if (y == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i= ini; i<= fin; i++)
    y->v[i - ini] = x->v[i];

  /* return data */
  return y;
}


/*** set all one ***/
void fxt_vecll_one(fxt_vecll *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_one: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = 1.0;
}


/*** element-wise negation ***/
void fxt_vecll_neg(fxt_vecll *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_neg: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = - x->v[i];
}


/*** element-wise square ***/
void fxt_vecll_sq(fxt_vecll *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_sq: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] *= x->v[i];
}


/*** element-wise square root ***/
void fxt_vecll_sqrt(fxt_vecll *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_sqrt: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    if (x->v[i] < 0.0) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecll_sqrt: negative entry\n");
      return;
    } else
      x->v[i] = fxt_lsqrt(x->v[i]);

  if (fxt_error_level() > 0)
    fxt_error_raise();
}


/*** scale a range ***/
void fxt_vecll_scale(fxt_vecll *x, long ini, long fin,
		     long double sc) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_scale: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_scale: irregal range\n");
    return;
  }

  /* scale the range */
  for (i= ini; i<= fin; i++)
    x->v[i] *= sc;
}


/*** element-wise multiply vector y to vector x ***/
void fxt_vecll_emul(fxt_vecll *x, fxt_vecll *y) {
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_emul: null pointer\n");
    return;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_emul: vectors different length\n");
    return;
  }

  /* multiplication */
  for (i=0; i< x->n; i++)
    x->v[i] *= y->v[i];
}


/*** inner product of two vectors ***/
long double fxt_vecll_dot(fxt_vecll *x, fxt_vecll *y) {
  long i;  long double s;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_dot: null pointer\n");
    return 0.0L;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_dot: size mismatch\n");
    return 0.0L;
  }

  /* multiplication */
  s = x->v[0] * y->v[0];

  for (i=1; i< x->n; i++)
    s += x->v[i] * y->v[i];

  /* return the value */
  return s;
}


/*** permute vector x according to perm ***/
void fxt_vecll_perm(fxt_vecll *x, fxt_vecl *perm) {
  long i, ii;
  long *idx;			/* original position */

  /* check null pointers */
  if (x == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_perm: null pointer\n");
    return;
  }

  /* check size */
  if (x->n != perm->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_perm: size mismatch\n");
    return;
  }

  /* check permutation vector */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecll_perm: error in permutation vector\n");
      return;
    }

  /* check very small cases */
  if (perm->n <= 1)
    return;

  /* allocate index array */
  idx = malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecll_perm: allocation failed\n");
    return;
  }

  /* check permutation vector again */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecll_perm: error in permutation vector\n");
      return;
    }

  /* initialize original index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute vector just as matrices */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the value does not come */
      /* save the current value */
      long double t = x->v[ii];

      i = ii;
      while (perm->v[i] != ii) {
	/* the index */
	long k = perm->v[i];

	/* move the value */
	x->v[i] = x->v[k];

	/* update index */
	idx[i] = k;

	/* goto next */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* move the value */
      x->v[i] = t;

      /* update index */
      idx[i] = ii;
    }

  /* deallocate temporal index array */
  free(idx);
}


/*** inverse permutation of x according to perm ***/
void fxt_vecll_iperm(fxt_vecll *x, fxt_vecl *perm) {
  long i, ii;
  long *idx;			/* original index */

  /* check null pointers */
  if (x == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_iperm: null pointer\n");
    return;
  }

  /* check size */
  if (x->n != perm->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_perm: vector size mismatch\n");
    return;
  }

  /* check permutation vector */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecll_perm: error in permutation vector\n");
      return;
    }

  /* check very small cases */
  if (perm->n <= 1)
    return;

  /* allocate index array */
  idx = malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecll_perm: allocation failed\n");
    return;
  }

  /* check permutation vector again */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecll_perm: error in permutation vector\n");
      return;
    }

  /* initialize original index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) {
      /* save the current value */
      long double t = x->v[ii];

      i = ii;
      while (perm->v[i] != ii) {
	/* the position */
	long k = perm->v[i];

	/* save the current value */
	long double tt = x->v[k];

	/* put the value */
	x->v[k] = t;

	/* update current value */
	t = tt;

	/* update index */
	idx[i] = k;

	/* goto next */
	i = k;
      }
      /* now perm->v[i] == ii */

      /* put the value */
      x->v[ii] = t;

      /* update index */
      idx[i] = ii;
    }

  /* deallocate temporal index array */
  free(idx);
}


/*** conversion from vecld ***/
fxt_vecll* fxt_vecld_to_vecll(fxt_vecld *x) {
  fxt_vecll *y;
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_to_vecll: null pointer\n");
    return NULL;
  }

  /* create new vector */
  y = fxt_vecll_new(x->n);
  if (y == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy the values */
  for (i=0; i< x->n; i++)
    y->v[i] = (long double) x->v[i];

  /* return result */
  return y;
}

/*** conversion to vecld ***/
fxt_vecld* fxt_vecll_to_vecld(fxt_vecll *x) {
  fxt_vecld *y;
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecll_to_vecld: null pointer\n");
    return NULL;
  }

  /* create new vector */
  y = fxt_vecld_new(x->n);
  if (y == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* copy values */
  for (i=0; i< x->n; i++)
    y->v[i] = (double) x->v[i];

  /* return result */
  return y;
}
