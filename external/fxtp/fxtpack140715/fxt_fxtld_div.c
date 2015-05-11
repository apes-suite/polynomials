#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  try division for FXT
*********************************************************************/

long fxtld_make_div(fxtld *fxt, fxt_matll *mat, fxt_vecll *x,
		    double eps, double th, long flop_c) {
  long nh, flop0, flop1;
  fxt_matll *hmat, *lmat;
  
  /* check null pointers */
  if (fxt == NULL || mat == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_make_div: null pointer\n");
    return 0;
  }

  /* negative upperbound */
  if (flop_c <= 0) {
    fxt->ch0 = fxt->ch1 = NULL;
    return 0;
  }

  /* half point */
  nh = (fxt->ui1 - fxt->ui0 + 1) / 2;

  /* make upper half */
  hmat = fxt_matll_clone(mat, 0, mat->nrow - 1, nh, mat->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* allocate data */
  fxt->ch1 = (fxtld*) malloc(sizeof(fxtld));
  if (fxt->ch1 == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxtld_make_div: allocation failed\n");
    return 0;
  }

  /* set parent */
  fxt->ch1->par = fxt;

  /* make ch1 */
  flop1 = fxtld_make(fxt->ch1, hmat, x, eps, th,
		     flop_c, fxt->ui0 + nh);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* check cost */
  if (flop1 >= flop_c) {
    /* flops over: undo make */
    fxtld_del(fxt->ch1);
    fxt->ch1 = NULL;
    fxt_matll_del(hmat);
    return flop1;
  }

  /* make lower half */
  lmat = fxt_matll_clone(mat, 0, mat->nrow - 1, 0, nh - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* allocate data */
  fxt->ch0 = (fxtld*) malloc(sizeof(fxtld));
  if (fxt->ch0 == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxtld_make_div: allocation failed\n");
    return 0;
  }

  /* set parent */
  fxt->ch0->par = fxt;

  /* make ch0 */
  flop0 = fxtld_make(fxt->ch0, lmat, x, eps, th,
		     flop_c - flop1, fxt->ui0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* check cost */
  if (flop0 + flop1 >= flop_c) {

    /* cost is over: undo both make */
    fxtld_del(fxt->ch0);
    fxt->ch0 = NULL;

    fxt_matll_del(lmat);

    fxtld_del(fxt->ch1);
    fxt->ch1 = NULL;

    fxt_matll_del(hmat);
  }

  return flop0 + flop1;
}
