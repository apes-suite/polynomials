#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  creation and deletion of FXTLD
*********************************************************************/

#define MAXITER 10

fxtld *fxtld_from_matll(fxt_matll *mat, fxt_vecll *x, double prec) {
  fxtld *fxt;			/* resulting fxt */
  double th, eps;
  long flop0, flop1;
  int c;

  /* check null pointers */
  if (mat == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_from_matll: null pointer\n");
    return NULL;
  }

  /* check size */
  if (mat->nrow != x->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_from_matll: size mismatch\n");
    return NULL;
  }

  /* safety factor with precision of the power method */
  prec *= 1.0 - POWERPREC;

  /* get the dropping threshould */
  th = fxt_matll_dropthreshould(mat, prec * 0.5);

  /* flops for direct computation */
  flop0 = fxt_matll_dropflop(mat, th);

  /* initial precision */
  eps = prec;

  for (c=0; c< MAXITER; c++) {
    double err, err2;
    fxt_matld *de;

    /* allocate memory */
    fxt = (fxtld*) malloc(sizeof(fxtld));
    if (fxt == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxtld_from_matll: allocation failed\n");
      return NULL;
    }

    /* the root structure */
    fxt->par = NULL;

    /* make fxt */
    flop1 = fxtld_make(fxt, mat, x, eps, th, flop0, 0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* get original matrix */
    de = fxt_matll_to_matld(mat);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* compute the error */
    err = fxtld_error(fxt, de, POWERPREC);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    /* for safety */
    err2 = fxtld_error(fxt, de, POWERPREC);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    if (err < err2)
      err = err2;

    /* deallocate temporal matrix */
    fxt_matld_del(de);

    /* check whether enough precision attained */
    if (err <= prec)
      break;

    /* delete data structure for retry */
    fxtld_del(fxt);
    fxt = NULL;

    /* decrease eps */
    if (err * 0.9 > prec)
      eps *= prec / err;
    else
      eps *= 0.9;
  }

  /* check success or not */
  if (fxt == NULL) {
#if DISALLOWHIGHERROR
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_from_matll: precision unattained\n");
#else
    fxt_error_set(FXT_ERROR_WARN,
		  "fxtld_from_matll: precision unattained\n");
#endif
    return NULL;
  }

  /* successful return */
  return fxt;
}


/*** deallocate data structrue ***/
void fxtld_del(fxtld *fxt) {
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_del: null pointer\n");
    return;
  }

  fxtld_unmake(fxt);
}

