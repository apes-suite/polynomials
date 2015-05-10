#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_vecld.h"

#include "fxt_vecl.h"
#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  vector of double, indexed with long int
*********************************************************************/

/*** create a new vector ***/
fxt_vecld* fxt_vecld_new(long size) {
  fxt_vecld *x;

  /* check input */
  if (size < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_new: negative size\n");
    return NULL;
  }

  /* allocate data structure */
  x = (fxt_vecld*) malloc(sizeof(fxt_vecld));
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecld_new: allocation failed\n");
    return NULL;
  }

  /* set size */
  x->n = size;

  /* allocate vector */
  if (size == 0)
    x->v = NULL;
  else {
    x->v = (double*) malloc(sizeof(double) * x->n);
    if (x->v == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_vecld_new: allocation failed\n");
      return NULL;
    }
  }

  /* return data */
  return x;
}


/*** deallocate memory ***/
void fxt_vecld_del(fxt_vecld *x) {

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_del: null pointer\n");
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
fxt_vecld* fxt_vecld_clone(fxt_vecld *x, long ini, long fin) {
  fxt_vecld *y;
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_clone: null pointer\n");
    return NULL;
  }

  /* check input */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_clone: irregal range\n");
    return NULL;
  }

  /* get new vector */
  y = fxt_vecld_new(fin - ini + 1);
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


/*** set part of vector ***/
void fxt_vecld_setprt(fxt_vecld *x, fxt_vecld *y, long ini, long fin){
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_setprt: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || y->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_setprt: irregal range\n");
    return;
  }

  /* check size */
  if (x->n != fin - ini + 1) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_setprt: size mismatch\n");
    return;
  }

  /* copy values */
  for (i=0; i< x->n; i++)
    x->v[i] = y->v[i + ini];
}


/*** set all zero ***/
void fxt_vecld_zero(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_zero: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = 0.0;
}


/*** set all one ***/
void fxt_vecld_one(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_one: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = 1.0;
}


/*** set random values ***/
void fxt_vecld_rand(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_rand: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = rand() / (1.0 + RAND_MAX);
}


/*** set k-th unit vector ***/
void fxt_vecld_unit(fxt_vecld *x, long k) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_unit: null pointer\n");
    return;
  }

  /* check input */
  if (k < 0 || x->n <= k) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_unit: irregal index\n");
    return;
  }

  /* set all zero */
  for (i=0; i< x->n; i++)
    x->v[i] = 0.0;

  /* set unit */
  x->v[k] = 1.0;
}


/*** element-wise negation ***/
void fxt_vecld_neg(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_neg: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] = - x->v[i];
}


/*** element-wise inversion ***/
void fxt_vecld_einv(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_einv: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    if (x->v[i] == 0.0) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_einv: divide by zero\n");
      return;
    } else
      x->v[i] = 1.0 / x->v[i];
}


/*** element-wise square ***/
void fxt_vecld_sq(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_sq: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    x->v[i] *= x->v[i];
}


/*** element-wise square root ***/
void fxt_vecld_sqrt(fxt_vecld *x) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_sqrt: null pointer\n");
    return;
  }

  for (i=0; i< x->n; i++)
    if (x->v[i] < 0.0) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_sqrt: negative entry\n");
      return;
    } else
      x->v[i] = sqrt(x->v[i]);
}


/*** scale a range ***/
void fxt_vecld_scale(fxt_vecld *x, long ini, long fin, double sc) {
  long i;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_scale: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_scale: irregal range\n");
    return;
  }

  /* scale the range */
  for (i= ini; i<= fin; i++)
    x->v[i] *= sc;
}


/*** add vector y to vector x ***/
void fxt_vecld_add(fxt_vecld *x, fxt_vecld *y) {
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_add: null pointer\n");
    return;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_add: size mismatch\n");
    return;
  }

  /* addition */
  for (i=0; i< x->n; i++)
    x->v[i] += y->v[i];
}


/*** subtract vector y from vector x ***/
void fxt_vecld_sub(fxt_vecld *x, fxt_vecld *y) {
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_sub: null pointer\n");
    return;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_sub: size mismatch\n");
    return;
  }

  /* subtraction */
  for (i=0; i< x->n; i++)
    x->v[i] -= y->v[i];
}


/*** element-wise multiply vector y to vector x ***/
void fxt_vecld_emul(fxt_vecld *x, fxt_vecld *y) {
  long i;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_emul: null pointer\n");
    return;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_emul: size mismatch\n");
    return;
  }

  /* multiplication */
  for (i=0; i< x->n; i++)
    x->v[i] *= y->v[i];
}


/*** inner product of two vectors ***/
double fxt_vecld_dot(fxt_vecld *x, fxt_vecld *y) {
  long i;  double s;

  /* check null pointers */
  if (x == NULL || y == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_dot: null pointer\n");
    return 0.0;
  }

  /* check input */
  if (x->n != y->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_dot: size mismatch\n");
    return 0.0;
  }

  /* multiplication */
  s = x->v[0] * y->v[0];

  for (i=1; i< x->n; i++)
    s += x->v[i] * y->v[i];

  /* return the value */
  return s;
}


/*** 2-norm ***/
double fxt_vecld_norm2(fxt_vecld *x) {

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_norm2: null pointer\n");
    return 0.0;
  }

  return sqrt(fxt_vecld_dot(x, x));
}


/*** normalize in 2-norm ***/
void fxt_vecld_normalize(fxt_vecld *x) {
  double s;

  /* check null pointers */
  if (x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_normalize: null pointer\n");
    return;
  }

  /* get the 2-norm */
  s = fxt_vecld_norm2(x);

  /* check zero vector */
  if (s == 0.0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_normalize: zero vector\n");
    return;
  }

  /* scale */
  fxt_vecld_scale(x, 0, x->n - 1, 1.0 / s);
}


/*** permute vector x according to perm ***/
void fxt_vecld_perm(fxt_vecld *x, fxt_vecl *perm) {
  long i, ii;
  long *idx;			/* original position */

  /* check null pointers */
  if (x == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_perm: null pointer\n");
    return;
  }

  /* check size */
  if (x->n != perm->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_perm: size mismatch\n");
    return;
  }

  /* check permutation vector */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_perm: error in permutation vector\n");
      return;
    }

  /* check very small cases */
  if (perm->n <= 1)
    return;

  /* allocate index array */
  idx = malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecld_perm: allocation failed\n");
    return;
  }

  /* check permutation vector again */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_perm: error in permutation vector\n");
      return;
    }

  /* initialize original index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  /* permute vector just as matrices */
  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) { /* the value does not come */
      /* save the current value */
      double t = x->v[ii];

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
void fxt_vecld_iperm(fxt_vecld *x, fxt_vecl *perm) {
  long i, ii;
  long *idx;			/* original index */

  /* check null pointers */
  if (x == NULL || perm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_iperm: null pointer\n");
    return;
  }

  /* check size */
  if (x->n != perm->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_perm: vector size mismatch\n");
    return;
  }

  /* check permutation vector */
  for (i=0; i< perm->n; i++)
    if (perm->v[i] < 0 || perm->n <= perm->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_perm: error in permutation vector\n");
      return;
    }

  /* check very small cases */
  if (perm->n <= 1)
    return;

  /* allocate index array */
  idx = malloc(sizeof(long) * perm->n);
  if (idx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_vecld_perm: allocation failed\n");
    return;
  }

  /* check permutation vector again */
  for (i=0; i< perm->n; i++)
    idx[i] = 0;

  for (i=0; i< perm->n; i++)
    if ((++ idx[perm->v[i]]) != 1) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_vecld_perm: error in permutation vector\n");
      return;
    }

  /* initialize original index */
  for (i=0; i< perm->n; i++)
    idx[i] = i;

  for (ii = 0; ii < perm->n; ii ++)
    if (idx[ii] != perm->v[ii]) {
      /* save the current value */
      double t = x->v[ii];

      i = ii;
      while (perm->v[i] != ii) {
	/* the position */
	long k = perm->v[i];

	/* save the current value */
	double tt = x->v[k];

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


/*** save a vector to a file ***/
void fxt_vecld_save(FILE *f, fxt_vecld *x) {

  /* check null pointers */
  if (f == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_save: null pointer\n");
    return;
  }

  fxt_file_writelongs(f, &x->n, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_file_writedoubles(f, x->v, x->n);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}


/*** restore a vector from a file ***/
fxt_vecld* fxt_vecld_restore(FILE *f) {
  long n;
  fxt_vecld *x;

  /* check null pointers */
  if (f == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_restore: null pointer\n");
    return NULL;
  }

  fxt_file_readlongs(f, &n, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (n < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_vecld_restore: negative size\n");
    return NULL;
  }

  x = fxt_vecld_new(n);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_file_readdoubles(f, x->v, n);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  return x;
}
