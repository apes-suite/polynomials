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

  dropping small entries
*********************************************************************/

/*** compute the threshold value for dropping ***/
double fxt_matld_dropthreshould(fxt_matld *mat, double eps) {
  double th;			/* the threshould */
  int c;			/* loop counter */

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_dropthreshold: null pointer\n");
    return 0.0;
  }

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
	toterr += v * v;

	/* update maximum error */
	if (maxerr < v)
	  maxerr = v;
      }
    }

    /* Frobenius norm */
    toterr = sqrt(toterr);

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


/*** flop counts when dropping ***/
long fxt_matld_dropflop(fxt_matld *mat, double th) {
  long flop = 0;		/* number of multiplications */
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_dropflop: null pointer\n");
    return 0;
  }

  /* count droppings */
  for (i=0; i < mat->nrow * mat->ncol; i++)
    if (fabs(mat->v[i]) > th)
      flop ++;

  return flop;
}


/*** do dropping ***/
void fxt_matld_drop(fxt_matld *mat, double th) {
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_drop: null pointer\n");
    return;
  }

  /* drop entries */
  for (i=0; i< mat->nrow * mat->ncol; i++)
    if (fabs(mat->v[i]) <= th)
      mat->v[i] = 0.0;
}
