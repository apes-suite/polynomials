#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fxt_error.h"
#include "fxt_sarld.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  creation and deletion of simple array representation
*********************************************************************/

static long ntmp, nent;

static void indexing(fxt_lagld_var *var, long nn) {
  fxt_lagld_ent *ent;

  if ((++ var->c) < var->n_out)
    return;

  var->c = 0;

  if (var->io == 0)
    var->v = (double) ((ntmp ++) + nn);

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    indexing(ent->ivar, nn), nent ++;
}
    

static void setarray(fxt_lagld_var *var, fxt_sarld *ar) {
  fxt_lagld_ent *ent;

  if ((++ var->c) < var->n_out)
    return;

  var->c = 0;

  if (var->io == 0) {
    long i = (long) var->v - ar->nrow;
    long j = ar->p[i];

    for (ent = var->oent; ent != NULL; ent = ent->onext) {
      ar->v[j] = ent->v;
      ar->k[j] = (long) ent->ovar->v;
      j ++;
    }

    ar->p[i + 1] = j;
  }

  for (ent = var->ient; ent != NULL; ent = ent->inext)
    setarray(ent->ivar, ar);
}


fxt_sarld* fxt_sarld_from_lagld(fxt_lagld* lag) {
  fxt_sarld *ar;  long i;

  for (i=0; i< lag->ncol; i++)
    if (lag->ivec->v[i].n_out == 0) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_sarld_new: column dropping\n");
      return NULL;
    }

  ntmp = nent = 0;
  for (i=0; i< lag->nrow; i++) {
    lag->ovec->v[i].v = (double) (lag->ncol + i);
    indexing(&lag->ovec->v[i], lag->ncol + lag->nrow);
  }

  ar = fxt_sarld_new(lag->nrow, lag->ncol, ntmp, nent);

  for (i=0; i< lag->ncol; i++) {
    fxt_lagld_var *var = &lag->ivec->v[i];
    fxt_lagld_ent *ent;
    long j = ar->p[i];

    for (ent = var->oent; ent != NULL; ent = ent->onext) {
      ar->v[j] = ent->v;
      ar->k[j] = (long) ent->ovar->v;
      j ++;
    }

    ar->p[i+1] = j;
  }
      
  for (i=0; i< lag->nrow; i++)
    setarray(&lag->ovec->v[i], ar);

  if (ar->p[ar->ncol + ar->ntmp] != nent) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_sarld_from_lagld: nent size mismatch\n");
    return NULL;
  }

  return ar;
}
