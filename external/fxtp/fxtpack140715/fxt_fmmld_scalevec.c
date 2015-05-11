#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  compute scaling vectors
*********************************************************************/


/*** compute scaling for multipole expansion ***/
static void comp_scalemul(fmmld_reg *reg, fmmld_reg *mreg,
			  fxt_matld *m0, fxt_vecld *v) {
  fxt_matld *m;

  if (mreg != NULL && (reg->type & FMM_RTMUL)) {
    m = fxt_matld_new(m0->nrow, reg->mm->ncol);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_matld_mul(m, m0, reg->mm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  } else
    m = m0;

  if (reg->type & FMM_RTMUL)
    mreg = reg;

  if (reg->type & FMM_RTEXP) {
    fxt_matld *m1;  fxt_vecld *v1;

    m1 = fxt_matld_new(m->nrow, reg->pm->ncol);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_matld_mul(m1, m, reg->pm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    v1 = fxt_vecld_new(m1->nrow);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_matld_norm2rows(v1, m1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_vecld_sq(v1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_vecld_add(v, v1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_vecld_del(v1);
    fxt_matld_del(m1);
  } else if (reg->ch0 != NULL) {
    comp_scalemul(reg->ch0, mreg, m, v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    comp_scalemul(reg->ch1, mreg, m, v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  if (m != m0)
    fxt_matld_del(m);
}


/*** get the scaling vector ***/
fxt_vecld* fxt_fmmld_get_scalemul(fmmld_reg *reg) {
  fxt_vecld *v;
  fxt_matld *m;

  /* check null pointers */
  if (reg == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_get_scalemul: null pointer\n");
    return NULL;
  }

  v = fxt_vecld_new(reg->km);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_zero(v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  m = fxt_matld_new(reg->km, reg->km);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_unit(m);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  comp_scalemul(reg, NULL, m, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_sqrt(v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_del(m);

  return v;
}



static void comp_scaleloc(fmmld_reg *reg,
			  fxt_matld *m0, fxt_vecld *v) {
  fmmld_cell *cell;

  if ((reg->type & FMM_RTLOC) && reg->ll != NULL) {
    fmmld_reg *lreg;  fxt_matld *m;

    for (lreg = reg->preg;
	 (lreg->type & FMM_RTLOC) == 0;
	 lreg = lreg->preg);

    m = fxt_matld_new(m0->nrow, reg->ll->ncol);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_matld_mul(m, m0, reg->ll);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    comp_scaleloc(lreg, m, v);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_matld_del(m);
  }

  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next)
    if (cell->type == FMM_CTTRN) {
      fxt_matld *m;

      m = fxt_matld_new(m0->nrow, cell->m->ncol);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_matld_mul(m, m0, cell->m);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      comp_scalemul(cell->evl_reg, NULL, m, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_matld_del(m);
    } else if (cell->type == FMM_CTLOC) {
      fxt_matld *m1;  fxt_vecld *v1;

      m1 = fxt_matld_new(m0->nrow, cell->m->ncol);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_matld_mul(m1, m0, cell->m);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      v1 = fxt_vecld_new(m1->nrow);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_matld_norm2rows(v1, m1);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_vecld_sq(v1);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_vecld_add(v, v1);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_vecld_del(v1);
      fxt_matld_del(m1);
    }
}


/*** get the scaling vector for local expansion ***/
fxt_vecld* fxt_fmmld_get_scaleloc(fmmld_reg *reg) {
  fxt_vecld *v;  fxt_matld *m;

  /* check null pointers */
  if (reg == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_get_scaleloc: null pointer\n");
    return NULL;
  }

  v = fxt_vecld_new(reg->kl);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_zero(v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  m = fxt_matld_new(reg->kl, reg->kl);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_unit(m);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  comp_scaleloc(reg, m, v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_sqrt(v);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_del(m);

  return v;
}
