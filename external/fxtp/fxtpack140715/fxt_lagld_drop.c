#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  compute weights
*********************************************************************/

static void weight_ivec(fxt_lagld_var *var, double a) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_in)
    return;

  for (ent = var->oent; ent != NULL; ent = ent->onext)
    weight_ivec(ent->ovar, ent->v * var->v);

  var->iw += var->v * var->v;
  var->c = 0;
}


static void weight_ovec(fxt_lagld_var *var, double a) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_out)
    return;

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    weight_ovec(ent->ivar, ent->v * var->v);

  var->ow += var->v * var->v;
  var->c = 0;
}


void fxt_lagld_setweight(fxt_lagld *lag) {
  long i, j;
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_setweight: null pointer\n");
    return;
  }

  for (i=0; i< lag->nu_var; i++)
    lag->valloc[i].iw = lag->valloc[i].ow = 0.0;

  for (i=0; i< lag->ncol; i++)
    for (j=0; j< lag->ncol; j++)
      weight_ivec(&lag->ivec->v[j], (i==j ? 1.0 : 0.0));

  for (i=0; i< lag->nrow; i++)
    for (j=0; j< lag->nrow; j++)
      weight_ovec(&lag->ovec->v[j], (i==j ? 1.0 : 0.0));

  for (i=0; i< lag->nu_ent; i++) {
    fxt_lagld_ent *ent = &lag->ealloc[i];
    if (ent->ivar != NULL && ent->ovar != NULL)
      ent->w = sqrt(ent->ivar->iw * ent->ovar->ow) * fabs(ent->v);
  }
}


static void add_var_evl(fxt_lagld_var *var, double a, double th) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_in)
    return;

  for (ent = var->oent; ent != NULL; ent = ent->onext)
    add_var_evl(ent->ovar,
		(ent->w <= th ? 0.0 : ent->v * var->v), th);

  var->c = 0;
}


void fxt_lagld_thevl(fxt_vecld *v, fxt_lagld *lag, fxt_vecld *u,
		     double th) {
  long i;
  
  /* check null pointers */
  if (v == NULL || lag == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_thevl: null pointer\n");
    return;
  }

  if (lag->ncol != u->n || lag->nrow != v->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_thevl: size mismatch\n");
    return;
  }

  for (i=0; i< lag->ncol; i++)
    add_var_evl(&lag->ivec->v[i], u->v[i], th);

  for (i=0; i< lag->nrow; i++)
    if (lag->ovec->v[i].n_in == 0)
      v->v[i] = 0.0;
    else
      v->v[i] = lag->ovec->v[i].v;
}


static void add_var_exp(fxt_lagld_var *var, double a, double th) {
  fxt_lagld_ent *ent;

  if (var->c ++ == 0)
    var->v = a;
  else
    var->v += a;

  if (var->c < var->n_out)
    return;

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    add_var_exp(ent->ivar,
		(ent->w <= th ? 0.0 : ent->v * var->v), th);

  var->c = 0;
}


void fxt_lagld_thexp(fxt_vecld *u, fxt_lagld *lag, fxt_vecld *v,
		     double th) {
  long i;
  
  /* check null pointers */
  if (u == NULL || lag == NULL || v == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_thexp: null pointer\n");
    return;
  }

  if (lag->ncol != u->n || lag->nrow != v->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_thexp: size mismatch\n");
    return;
  }

  for (i=0; i< lag->nrow; i++)
    add_var_exp(&lag->ovec->v[i], v->v[i], th);

  for (i=0; i< lag->ncol; i++)
    if (lag->ivec->v[i].n_out == 0)
      u->v[i] = 0.0;
    else
      u->v[i] = lag->ivec->v[i].v;
}


void fxt_lagld_thzero(fxt_lagld *lag, double th) {
  long i;
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_thzero: null pointer\n");
    return;
  }

  for (i=0; i< lag->nu_ent; i++) {
    fxt_lagld_ent *ent = &lag->ealloc[i];

    if (ent->w <= th)
      ent->v = 0.0;
  }
}
