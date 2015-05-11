#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matll.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

   drop small entries
*********************************************************************/


/*** compute threshould for dropping ***/
double fxt_matll_dropthreshould(fxt_matll *mat, double eps) {
  double th;			/* threshould */
  fxt_matld *m;			/* double-version clone */

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_dropthreshould: null pointer\n");
    return 0.0;
  }

  /* make double-version clone */
  m = fxt_matll_to_matld(mat);
  if (m == NULL) {
    fxt_error_raise();
    return 0.0;
  }

  /* compute the threshould */
  th = fxt_matld_dropthreshould(m, eps);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  /* deallocate the clone */
  fxt_matld_del(m);

  /* return results */
  return th;
}


/*** the number of non-zero entries after dropping ***/
long fxt_matll_dropflop(fxt_matll *mat, double th) {
  long flop = 0;		/* number of multiplications */
  long i;

  /* check null pointers */
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_dropflop: null pointer\n");
    return 0;
  }

  /* count entries larger than th */
  for (i=0; i< mat->nrow * mat->ncol; i++)
    if (fabs((double) mat->v[i]) > th)
      flop ++;

  return flop;
}
