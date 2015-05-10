#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for matrix-vector product of double
*********************************************************************/

/*** general add matrix-vector product ***/
static void addmvprod(fxt_vecld *y, long yini, long yfin,
		      fxt_matld *a,
		      fxt_vecld *x, long xini, long xfin) {
  long i, j;

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_vecld: null pointer\n");
    return;
  }

  /* check range */
  if (0 > yini || yfin >= y->n || 0 > xini || xfin >= x->n) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: irregal range\n");
    return;
  }

  /* check size */
  if (a->nrow != yfin - yini + 1 || a->ncol != xfin - xini + 1) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: size mismatch\n");
    return;
  }

  /* compute the product */
  for (i=0; i< a->nrow; i++) {
    /* clear (ncol may be zero) */
    double s = 0.0;

    for (j=0; j< a->ncol; j++)
      s += MENT(a, i, j) * x->v[j + xini];

    y->v[i + yini] += s;
  }
}


/*** general set matrix-vector product ***/
static void setmvprod(fxt_vecld *y, long yini, long yfin,
		      fxt_matld *a,
		      fxt_vecld *x, long xini, long xfin) {
  long i;

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_vecld: null pointer\n");
    return;
  }

  /* check range */
  if (0 > yini || yfin >= y->n) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: irregal range\n");
    return;
  }

  /* clear y */
  for (i = yini; i <= yfin; i++)
    y->v[i] = 0.0;

  /* compute the product */
  addmvprod(y, yini, yfin, a, x, xini, xfin);
}


/*** set matrix-vector product y = a * x ***/
void fxt_matld_setmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setmatvec: null pointer\n");
    return;
  }

  /* check size */
  if (y->n != a->nrow || x->n != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  setmvprod(y, 0, y->n - 1, a, x, 0, x->n - 1);
}


/*** add matrix-vector product y += a * x ***/
void fxt_matld_addmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_admatvec: null pointer\n");
    return;
  }

  /* check size */
  if (y->n != a->nrow || x->n != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addmvprod(y, 0, y->n - 1, a, x, 0, x->n - 1);
}


/*** set matrix-partial vector product y = a * x[ini:fin] ***/
void fxt_matld_setmatprtvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x,
			    long ini, long fin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (y->n != a->nrow || fin - ini + 1 != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_setmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  setmvprod(y, 0, y->n - 1, a, x, ini, fin);
}


/*** add matrix-partial vector product y += a * x[ini:fin] ***/
void fxt_matld_addmatprtvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x,
			    long ini, long fin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (y->n != a->nrow || fin - ini + 1 != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addmatprtvec: size mismatch\n");
    return;
  }

  /* compute product */
  addmvprod(y, 0, y->n - 1, a, x, ini, fin);
}


/*** partial set matrix-vector product y[ini:fin] = a * x ***/
void fxt_matld_prtsetmatvec(fxt_vecld *y, long ini, long fin,
			    fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || y->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatvec: irregal range\n");
    return;
  }

  /* check size */
  if (fin - ini + 1 != a->nrow || x->n != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatvec: size mismatch\n");
    return;
  }

  /* compute product */
  setmvprod(y, ini, fin, a, x, 0, x->n - 1);
}


/*** partial add matrix-vector product y[ini:fin] += a * x ***/
void fxt_matld_prtaddmatvec(fxt_vecld *y, long ini, long fin,
			    fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || y->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatvec: irregal range\n");
    return;
  }

  /* check size */
  if (fin - ini + 1 != a->nrow || x->n != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatvec: size mismatc\n");
    return;
  }

  /* compute product */
  addmvprod(y, ini, fin, a, x, 0, x->n - 1);
}


/*** partial set matrix-partial vector product ***/
void fxt_matld_prtsetmatprtvec(fxt_vecld *y, long yini, long yfin,
			       fxt_matld *a,
			       fxt_vecld *x, long xini, long xfin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (yfin < yini || yini < 0 || y->n <= yfin ||
      xfin < xini || xini < 0 || x->n <= xfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (yfin - yini + 1 != a->nrow || xfin - xini + 1 != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsetmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  setmvprod(y, yini, yfin, a, x, xini, xfin);
}


/*** partial add matrix-partial vector product ***/
void fxt_matld_prtaddmatprtvec(fxt_vecld *y, long yini, long yfin,
			       fxt_matld *a,
			       fxt_vecld *x, long xini, long xfin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (yini < 0 || y->n <= yfin || xini < 0 || x->n <= xfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (yfin - yini + 1 != a->nrow || xfin - xini + 1 != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addmvprod(y, yini, yfin, a, x, xini, xfin);
}


/*** general add transposed matrix-vector product ***/
static void addtrmvprod(fxt_vecld *y, long yini, long yfin,
			fxt_matld *a,
			fxt_vecld *x, long xini, long xfin) {
  long i, j;

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_vecld: null pointer\n");
    return;
  }

  /* check range */
  if (0 > yini || yfin >= y->n || 0 > xini || xfin >= x->n) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: irregal range\n");
    return;
  }

  /* check size */
  if (yfin - yini + 1 != a->ncol || xfin - xini + 1 != a->nrow) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: size mismatch\n");
    return;
  }

  /* compute the product */
  for (i=0; i< a->nrow; i++)
    for (j=0; j< a->ncol; j++)
      y->v[j + yini] += MENT(a, i, j) * x->v[i + xini];
}


/*** general set transposed matrix-vector product ***/
static void settrmvprod(fxt_vecld *y, long yini, long yfin,
			fxt_matld *a,
			fxt_vecld *x, long xini, long xfin) {
  long i;

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_vecld: null pointer\n");
    return;
  }

  /* check range */
  if (0 > yini || yfin >= y->n) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_vecld: irregal range\n");
    return;
  }

  /* clear y */
  for (i= yini; i <= yfin; i++)
    y->v[i] = 0.0;

  /* compute the product */
  addtrmvprod(y, yini, yfin, a, x, xini, xfin);
}


/*** set transposed matrix-vector product y = A^t * x ***/
void fxt_matld_settrmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_settrmatvec: null pointer\n");
    return;
  }

  /* check size */
  if (y->n != a->ncol || x->n != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_settrmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  settrmvprod(y, 0, y->n - 1, a, x, 0, x->n - 1);
}


/*** add transposed matrix-vector product y += A^t * x ***/
void fxt_matld_addtrmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addtrmatvec: null pointer\n");
    return;
  }

  /* check size */
  if (y->n != a->ncol || x->n != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addtrmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addtrmvprod(y, 0, y->n - 1, a, x, 0, x->n - 1);
}


/*** y = A^t * x[ini:fin] ***/
void fxt_matld_settrmatprtvec(fxt_vecld *y, fxt_matld *a,
			      fxt_vecld *x, long ini, long fin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_settrmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_settrmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (y->n != a->ncol || fin - ini + 1 != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_settrmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  settrmvprod(y, 0, y->n - 1, a, x, ini, fin);
}


/*** y += A^t * x[ini:fin] ***/
void fxt_matld_addtrmatprtvec(fxt_vecld *y, fxt_matld *a,
			      fxt_vecld *x, long ini, long fin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addtrmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || x->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addtrmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (y->n != a->ncol || fin - ini + 1 != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_addtrmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addtrmvprod(y, 0, y->n - 1, a, x, ini, fin);
}


/*** y[ini:fin] = A^t * x ***/
void fxt_matld_prtsettrmatvec(fxt_vecld *y, long ini, long fin,
			      fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || y->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatvec: irregal range\n");
    return;
  }

  /* check size */
  if (fin - ini + 1 != a->ncol || x->n != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  settrmvprod(y, ini, fin, a, x, 0, x->n - 1);
}


/*** y[ini:fin] += A^t * x ***/
void fxt_matld_prtaddtrmatvec(fxt_vecld *y, long ini, long fin,
			      fxt_matld *a, fxt_vecld *x) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatvec: null pointer\n");
    return;
  }

  /* check range */
  if (fin < ini || ini < 0 || y->n <= fin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatvec: irregal range\n");
    return;
  }

  /* check size */
  if (fin - ini + 1 != a->ncol || x->n != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addtrmvprod(y, ini, fin, a, x, 0, x->n - 1);
}


/*** y[yini:yfin] = A^t * x[xini:xfin] ***/
void fxt_matld_prtsettrmatprtvec(fxt_vecld *y, long yini, long yfin,
				 fxt_matld *a,
				 fxt_vecld *x, long xini, long xfin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (yfin < yini || yini < 0 || y->n <= yfin ||
      xfin < xini || xini < 0 || x->n <= xfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (yfin - yini + 1 != a->ncol || xfin - xini + 1 != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtsettrmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  settrmvprod(y, yini, yfin, a, x, xini, xfin);
}


/*** y[yini:yfin] += A^t * x[xini:xfin] ***/
void fxt_matld_prtaddtrmatprtvec(fxt_vecld *y, long yini, long yfin,
				 fxt_matld *a,
				 fxt_vecld *x, long xini, long xfin) {

  /* check null pointers */
  if (y == NULL || a == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatprtvec: null pointer\n");
    return;
  }

  /* check range */
  if (yfin < yini || yini < 0 || y->n <= yfin ||
      xfin < xini || xini < 0 || x->n <= xfin) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatprtvec: irregal range\n");
    return;
  }

  /* check size */
  if (yfin - yini + 1 != a->ncol || xfin - xini + 1 != a->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_prtaddtrmatprtvec: size mismatch\n");
    return;
  }

  /* compute the product */
  addtrmvprod(y, yini, yfin, a, x, xini, xfin);
}
