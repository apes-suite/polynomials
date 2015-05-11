#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  linear algorithm graph
*********************************************************************/

fxt_lagld* fxt_lagld_new(long ncol, long nrow, long nvar, long nent) {
  fxt_lagld *lag;
  long i;

  /* check input */
  if (ncol < 0 || nrow < 0 || nvar < 0 || nent < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_new: negative size\n");
    return NULL;
  }

  /* allocate memory */
  lag = (fxt_lagld*) malloc(sizeof(fxt_lagld));
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_lagld_new: allocation failed\n");
    return NULL;
  }

  /* initialize entries */
  lag->ncol = ncol;
  lag->nrow = nrow;

  lag->nu_var = 0;
  lag->nu_ent = 0;

  /* allocate variables */
  lag->na_var = ncol + nrow + nvar;

  lag->valloc = (fxt_lagld_var*)
    malloc(sizeof(fxt_lagld_var) * lag->na_var);
  if (lag->valloc == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_lagld_new: allocation failed\n");
    return NULL;
  }

  /* initialize variables */
  for (i=0; i< lag->na_var; i++) {
    fxt_lagld_var *var = &lag->valloc[i];

    /* default to internal variable */
    var->io = 0;

    /* number of in/out */
    var->n_in = var->n_out = 0;

    /* lists */
    var->ient = var->oent = NULL;

    /* value and counter */
    var->v = 0.0;
    var->c = 0;

    /* weights */
    var->iw = var->ow = 0.0;
  }

  /* allocate entries */
  lag->na_ent = nent;

  lag->ealloc = (fxt_lagld_ent*)
    malloc(sizeof(fxt_lagld_ent) * lag->na_ent);
  if (lag->ealloc == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_lagld_new: allocation failed\n");
    return NULL;
  }

  /* make input vector */
  lag->ivec = fxt_lagld_newvec(lag, ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* mark as input vector */
  for (i=0; i< ncol; i++) {
    lag->ivec->v[i].io = 1;
    lag->ivec->v[i].n_in = 1;
  }

  /* make output vector */
  lag->ovec = fxt_lagld_newvec(lag, nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* mark as output vector */
  for (i=0; i< nrow; i++) {
    lag->ovec->v[i].io = 2;
    lag->ovec->v[i].n_out = 1;
  }

  return lag;
}


void fxt_lagld_del(fxt_lagld *lag) {
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_del: null pointer\n");
    return;
  }

  free(lag->ovec);
  lag->ovec = NULL;

  free(lag->ivec);
  lag->ivec = NULL;

  free(lag->ealloc);
  lag->ealloc = NULL;

  free(lag->valloc);
  lag->valloc = NULL;

  free(lag);
}
