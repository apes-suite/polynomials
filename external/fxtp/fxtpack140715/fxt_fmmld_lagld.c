#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fxt_error.h"
#include "fxt_fmmld_loc.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  FMM to LAG data structure
*********************************************************************/


static void prep_lagld(fxt_lagld *lag, fmmld_reg *reg) {
  if ((reg->type & FMM_RTMUL) && reg->km > 0) {

    reg->vm = fxt_lagld_newvec(lag, reg->km);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else
    reg->vm = NULL;

  if ((reg->type & FMM_RTLOC) && reg->kl > 0) {

    reg->vl = fxt_lagld_newvec(lag, reg->kl);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else
    reg->vl = NULL;

  if (reg->ch0 != NULL) {

    prep_lagld(lag, reg->ch0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    prep_lagld(lag, reg->ch1);
  }
}


static void mul_lagld(fxt_lagld *lag, fmmld_reg *reg,
		      fxt_lagld_vec *u, fmmld_reg *preg) {
  fmmld_reg *mreg;

  if (reg->type & FMM_RTMUL)
    mreg = reg;
  else
    mreg = preg;

  if ((reg->type & FMM_RTEXP) && mreg != NULL) {

    if (mreg->km > 0 && reg->ns > 0) {
      fxt_matld_lagld(lag, mreg->vm, 0, reg->pm, u, reg->is);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    mreg = NULL;
  }

  if (reg->ch0 != NULL) {

    mul_lagld(lag, reg->ch0, u, mreg);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    mul_lagld(lag, reg->ch1, u, mreg);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  }

  if ((reg->type & FMM_RTMUL) && preg != NULL) {

    if (preg->km > 0 && reg->km > 0) {
      fxt_matld_lagld(lag, preg->vm, 0, reg->mm, reg->vm, 0);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }
  }
}


static void loc_lagld(fxt_lagld *lag, fxt_lagld_vec *v,
		      fmmld_reg *reg, fxt_lagld_vec *u,
		      fmmld_reg *preg) {
  fmmld_cell *cell;

  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;

    if (cell->type == FMM_CTTRN) {

      if (reg->kl > 0 && creg->km > 0) {
	fxt_matld_lagld(lag, reg->vl, 0, cell->m, creg->vm, 0);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

    } else if (cell->type == FMM_CTLOC) {

      if (reg->kl > 0 && creg->ns > 0) {
	fxt_matld_lagld(lag, reg->vl, 0, cell->m, u, creg->is);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

    } else if (cell->type == FMM_CTDIR) {

      if (reg->nt > 0 && creg->ns > 0) {
	fxt_matld_lagld(lag, v, reg->it, cell->m, u, creg->is);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

    } else if (cell->type == FMM_CTMUL) {

      if (reg->nt > 0 && creg->km > 0) {
	fxt_matld_lagld(lag, v, reg->it, cell->m, creg->vm, 0);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

    }
  }

  if (reg->type & FMM_RTLOC) {

    if (preg != NULL) {

      if (reg->kl > 0 && preg->kl > 0) {
	fxt_matld_lagld(lag, reg->vl, 0, reg->ll, preg->vl, 0);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

    }

    preg = reg;
  }

  if ((reg->type & FMM_RTEVL) && preg != NULL) {

    if (reg->nt > 0 && preg->kl > 0) {
      fxt_matld_lagld(lag, v, reg->it, reg->lp, preg->vl, 0);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    preg = NULL;
  }

  if (reg->ch0 != NULL) {

    loc_lagld(lag, v, reg->ch0, u, preg);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    loc_lagld(lag, v, reg->ch1, u, preg);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  }
}


static void finl_lagld(fmmld_reg *reg) {

  if (reg->ch0 != NULL) {
    finl_lagld(reg->ch1);
    finl_lagld(reg->ch0);
  }

  if (reg->vl != NULL) {
    free(reg->vl);
    reg->vl = NULL;
  }

  if (reg->vm != NULL) {
    free(reg->vm);
    reg->vm = NULL;
  }
}


void fxt_fmmld_lagld(fxt_lagld *lag, fxt_lagld_vec *v,
		     fxt_fmmld *fmm, fxt_lagld_vec *u) {
  int i;

  /* check null pointers */
  if (lag == NULL || v == NULL || fmm == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_lagld: null pointer\n");
    return;
  }

  /* check size */
  if (v->n != fmm->nt || u->n != fmm->ns) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_lagld: size mismatch\n");
    return;
  }

  for (i=0; i< fmm->n_roots; i++) {
    prep_lagld(lag, &fmm->roots[i]);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  for (i=0; i< fmm->n_roots; i++) {
    mul_lagld(lag, &fmm->roots[i], u, NULL);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  for (i=0; i< fmm->n_roots; i++) {
    loc_lagld(lag, v, &fmm->roots[i], u, NULL);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  for (i = fmm->n_roots - 1; i >= 0; i--)
    finl_lagld(&fmm->roots[i]);
}
