#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  normalize linear algorithm graph
*********************************************************************/

static void del_ent(fxt_lagld_ent *ent) {
  fxt_lagld_var *var;

  /* remove from input variable */
  var = ent->ivar;

  var->n_out --;

  if (var->oent == ent)
    var->oent = ent->onext;

  else {
    fxt_lagld_ent *ee;

    for (ee = var->oent; ee->onext != ent; ee = ee->onext);

    ee->onext = ent->onext;
  }

  /* remove from output variable */
  var = ent->ovar;

  var->n_in --;

  if (var->ient == ent)
    var->ient = ent->inext;

  else {
    fxt_lagld_ent *ee;

    for (ee = var->ient; ee->inext != ent; ee = ee->inext);

    ee->inext = ent->inext;
  }

  ent->ivar = ent->ovar = NULL;
}


static void touch_forward(fxt_lagld_var *var) {
  fxt_lagld_ent *ent;

  if (var->c != 0)
    return;

  var->c = 1;

  for (ent = var->oent; ent != NULL; ent = ent->onext)
    if (ent->v != 0.0)
      touch_forward(ent->ovar);
}


static void touch_backward(fxt_lagld_var *var) {
  fxt_lagld_ent *ent;

  if (var->c != 1)
    return;

  var->c = 2;

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    if (ent->v != 0.0)
      touch_backward(ent->ivar);
}


static void remove_onein(fxt_lagld_var *var) {
  fxt_lagld_ent *ient, *oent, *lent;
  fxt_lagld_var *ivar;
  double v;

  if (var->n_in != 1 || var->io != 0)
    return;

  /* input entry: to be removed */
  ient = var->ient;

  /* the previous node */
  ivar = ient->ivar;

  /* the value of removed entry */
  v = ient->v;

  /* remove entry */
  del_ent(ient);

  /* output entry: to be attached to ivar */
  oent = var->oent;

  /* search for the last entry, scaling by v */
  for (lent = oent; lent->onext != NULL; lent = lent->onext) {
    lent->v *= v;
    lent->ivar = ivar;
  }

  lent->v *= v;
  lent->ivar = ivar;

  lent->onext = ivar->oent;
  ivar->oent = oent;

  ivar->n_out += var->n_out;

  var->oent = NULL;
  var->n_out = 0;
}


static void remove_oneout(fxt_lagld_var *var) {
  fxt_lagld_ent *oent, *ient, *lent;
  fxt_lagld_var *ovar;
  double v;

  if (var->n_out != 1 || var->io != 0)
    return;

  /* output entry: to be removed */
  oent = var->oent;

  /* next variable */
  ovar = oent->ovar;

  /* value of the entry */
  v = oent->v;

  /* remove entry */
  del_ent(oent);

  /* input list: to be attached to ovar */
  ient = var->ient;

  /* search for the last entry */
  for (lent = ient; lent->inext != NULL; lent = lent->inext) {
    lent->v *= v;
    lent->ovar = ovar;
  }

  lent->v *= v;
  lent->ovar = ovar;

  lent->inext = ovar->ient;
  ovar->ient = ient;

  ovar->n_in += var->n_in;

  var->ient = NULL;
  var->n_in = 0;
}


static void remove_twobytwo(fxt_lagld_var *var) {
  fxt_lagld_ent *enta, *entb, *entc, *entd;
  double va, vb, vc, vd;

  if (var->io != 0 || var->n_in != 2 || var->n_out != 2)
    return;
  
  enta = var->oent;
  entb = enta->onext;
  entc = var->ient;
  entd = entc->inext;

  va = enta->v;
  vb = entb->v;
  vc = entc->v;
  vd = entd->v;

  enta->v *= vc;		/* XA -> CA */
  entb->v *= vd;		/* XB -> DB */
  entc->v *= vb;		/* CX -> CB */
  entd->v *= va;		/* DX -> DA */

  enta->ivar = entc->ivar;
  enta->onext = enta->ivar->oent;
  enta->ivar->oent = enta;
  enta->ivar->n_out ++;

  entb->ivar = entd->ivar;
  entb->onext = entb->ivar->oent;
  entb->ivar->oent = entb;
  entb->ivar->n_out ++;

  entc->ovar = entb->ovar;
  entc->inext = entc->ovar->ient;
  entc->ovar->ient = entc;
  entc->ovar->n_in ++;

  entd->ovar = enta->ovar;
  entd->inext = entd->ovar->ient;
  entd->ovar->ient = entd;
  entd->ovar->n_in ++;

  var->ient = var->oent = NULL;
  var->n_in = var->n_out = 0;
}


void fxt_lagld_normalize(fxt_lagld *lag) {
  long i;
  
  /* check null pointers */
  if (lag == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_lagld_normalize: null pointer\n");
    return;
  }

  /* check reachability from inputs */
  for (i=0; i< lag->ncol; i++)
    touch_forward(&lag->ivec->v[i]);

  /* check reachability from outputs */
  for (i=0; i< lag->nrow; i++)
    touch_backward(&lag->ovec->v[i]);

  /* remove entries */
  for (i=0; i< lag->nu_ent; i++) {
    fxt_lagld_ent *ent = &lag->ealloc[i];

    if (ent->ivar == NULL || ent->ovar == NULL)
      continue;

    if (ent->v == 0.0 || ent->ivar->c != 2 || ent->ovar->c != 2)
      del_ent(ent);
  }

  /* clear mark */
  for (i=0; i< lag->nu_var; i++) {
    fxt_lagld_var *var = &lag->valloc[i];

    if (var->c != 2 && var->io == 0) {
      var->ient = var->oent = NULL;
      var->n_in = var->n_out = 0;
    }

    var->c = 0;
  }

  /* remove oneout */
  for (i=0; i< lag->nu_var; i++) {
    fxt_lagld_var *var = &lag->valloc[i];

    if (var->io == 0 && var->n_out == 1)
      remove_oneout(var);
  }

  /* remove onein */
  for (i=0; i< lag->nu_var; i++) {
    fxt_lagld_var *var = &lag->valloc[i];

    if (var->io == 0 && var->n_in == 1)
	remove_onein(var);
  }

  /* remove twobytwo */
  for (i=0; i< lag->nu_var; i++) {
    fxt_lagld_var *var = &lag->valloc[i];

    if (var->io == 0 && var->n_in == 2 && var->n_out == 2)
      remove_twobytwo(var);
  }
}
