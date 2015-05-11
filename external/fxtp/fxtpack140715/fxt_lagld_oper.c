#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  basic operations for linear algorithm graph
*********************************************************************/

fxt_lagld_vec* fxt_lagld_ivec(fxt_lagld *lag) {
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_ivec: null pointer\n");
    return NULL;
  }

  return lag->ivec;
}


fxt_lagld_vec* fxt_lagld_ovec(fxt_lagld *lag) {
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_ovec: null pointer\n");
    return NULL;
  }

  return lag->ovec;
}


fxt_lagld_vec* fxt_lagld_newvec(fxt_lagld *lag, long n) {
  fxt_lagld_vec *v;
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_newvec: null pointer\n");
    return NULL;
  }

  /* check size */
  if (n < 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_newvec: negative size\n");
    return NULL;
  }

  /* check allocation */
  if (lag->nu_var + n > lag->na_var) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_newvec: allocation over\n");
    return NULL;
  }

  /* allocate structure */
  v = (fxt_lagld_vec*) malloc(sizeof(fxt_lagld_vec));
  if (v == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_lagld_newvec: allocation failed\n");
    return NULL;
  }

  /* set fields */
  v->n = n;

  v->v = &lag->valloc[lag->nu_var];
  lag->nu_var += n;

  return v;
}


void fxt_lagld_add(fxt_lagld *lag, fxt_lagld_var *y,
		   double a, fxt_lagld_var *x) {
  fxt_lagld_ent *e;
  
  /* check null pointers */
  if (lag == NULL || y == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_add: null pointer\n");
    return;
  }

  /* search for the same pair */
  for (e = y->ient; e != NULL && e->ivar != y; e = e->inext);

  if (e != NULL) {
    /* same pair found */

    e->v += a;
    return;

  }
 
  /* check allocated size */
  if (lag->nu_ent >= lag->na_ent) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_add: allocation over\n");
    return;
  }

  /* allocate new entry */
  e = &lag->ealloc[lag->nu_ent ++];

  /* set value */
  e->v = a;

  /* set output variable */
  e->ovar = y;
  y->n_in ++;
  e->inext = y->ient;
  y->ient = e;

  /* set input variable */
  e->ivar = x;
  x->n_out ++;
  e->onext = x->oent;
  x->oent = e;

  /* initialize */
  e->w = 0.0;
}
