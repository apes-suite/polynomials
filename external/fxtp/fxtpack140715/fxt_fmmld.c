#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

#include "fxt_vecld.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generalized one-dimensional fast multipole method
*********************************************************************/

/*** create a new object ***/
fxt_fmmld* fxt_fmmld_new(fxt_matld *mat, fxt_vecld *xs, fxt_vecld *xt,
			 double prec) {
  fxt_fmmld *fmm;
  long i;

  /* check null pointers */
  if (mat == NULL || xs == NULL || xt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_new: null pointer\n");
    return NULL;
  }

  /* check size */
  if (mat->ncol != xs->n || mat->nrow != xt->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_new: size mismatch\n");
    return NULL;
  }

  /* check xs */
  for (i=1; i< xs->n; i++)
    if (xs->v[i-1] >= xs->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_fmmld_new: unsorted source points\n");
      return NULL;
    }

  /* check xt */
  for (i=1; i< xt->n; i++)
    if (xt->v[i-1] >= xt->v[i]) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_fmmld_new: unsorted target points\n");
      return NULL;
    }

  /* allocate data structure */
  fmm = (fxt_fmmld*) malloc(sizeof(fxt_fmmld));
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_fmmld_new: allocation failed\n");
    return NULL;
  }

  /* definition parameters */
  fmm->x = xs;
  fmm->y = xt;
  fmm->mat = mat;

  fmm->ns = mat->ncol;
  fmm->nt = mat->nrow;
  fmm->prec = prec;

  fmm->cnorm = fxt_matld_norm2(mat, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* clear data */
  fmm->roots = NULL;

  fmm->kest = fmm->kmax = 0;
  fmm->min_np = 0;

  fmm->warn = 0;

  /* make FMM */
  fxt_fmmld_preproc(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* return result */
  return fmm;
}


/*** returns non-zero if the precision unattained ***/
int fxt_fmmld_warn(fxt_fmmld *fmm) {

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_warn: null pointer\n");
    return 0;
  }

  return fmm->warn;
}


/*** deallocate memory ***/
void fxt_fmmld_del(fxt_fmmld *fmm) {

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_del: null pointer\n");
    return;
  }

  /* undo preprocessing */
  fxt_fmmld_unpreproc(fmm);

  /* unlink definitions for safety */
  fmm->x = fmm->y = NULL;
  fmm->mat = NULL;

  /* deallocate body */
  free(fmm);
}
