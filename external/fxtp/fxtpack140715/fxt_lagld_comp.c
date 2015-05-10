#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  evaluate linear algorithm graph
*********************************************************************/

static void add_var_evl(fxt_lagld_var *var, double a) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_in)
    return;

  /* additions complete */

  var->c = 0;

  for (ent = var->oent; ent != NULL; ent = ent->onext)
    add_var_evl(ent->ovar, ent->v * var->v);
}


void fxt_lagld_evl(fxt_vecld *v, fxt_lagld *lag, fxt_vecld *u) {
  long i;
  
  /* check null pointers */
  if (v == NULL || lag == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_evl: null pointer\n");
    return;
  }

  /* check size */
  if (v->n != lag->nrow || u->n != lag->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_evl: size mismatch\n");
    return;
  }

  /* set input */
  for (i=0; i< lag->ncol; i++)
    add_var_evl(&lag->ivec->v[i], u->v[i]);

  /* get output */
  for (i=0; i< lag->nrow; i++)
    if (lag->ovec->v[i].n_in == 0)
      v->v[i] = 0.0;
    else
      v->v[i] = lag->ovec->v[i].v;
}


static void add_var_exp(fxt_lagld_var *var, double a) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_out)
    return;

  /* addition complete */

  var->c = 0;

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    add_var_exp(ent->ivar, ent->v * var->v);
}


void fxt_lagld_exp(fxt_vecld *u, fxt_lagld *lag, fxt_vecld *v) {
  long i;
  
  /* check null pointers */
  if (u == NULL || lag == NULL || v == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_exp: null pointer\n");
    return;
  }

  /* check size */
  if (u->n != lag->ncol || v->n != lag->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_exp: size mismatch\n");
    return;
  }

  /* set input */
  for (i=0; i< lag->nrow; i++)
    add_var_exp(&lag->ovec->v[i], v->v[i]);

  /* get output */
  for (i=0; i< lag->ncol; i++)
    if (lag->ivec->v[i].n_out == 0)
      u->v[i] = 0.0;
    else
      u->v[i] = lag->ivec->v[i].v;
}


static long mflop;


static void count_mflop(fxt_lagld_var *var) {
  fxt_lagld_ent *ent;

  if (++ var->c < var->n_in)
    return;

  var->c = 0;

  for (ent = var->oent; ent != NULL; ent = ent->onext) {
    mflop ++;
    count_mflop(ent->ovar);
  }
}


long fxt_lagld_evaluate_mflop(fxt_lagld *lag) {
  long i;
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_evaluate_mflop: null pointer\n");
    return 0;
  }

  mflop = 0;

  for (i=0; i< lag->ncol; i++)
    count_mflop(&lag->ivec->v[i]);

  return mflop;
}
