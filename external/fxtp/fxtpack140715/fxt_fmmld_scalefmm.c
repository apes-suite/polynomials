#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  scale columns and rows of the FMM
*********************************************************************/

/*** source scaling ***/
static void scalecol_reg(fmmld_reg *reg, fxt_vecld *s) {
  fmmld_cell *cell;
  fxt_vecld *ss;

  if (reg->ns == 0)
    return;

  /* the corresponding part of scaling vector */
  ss = fxt_vecld_clone(s, reg->is, reg->is + reg->ns - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* scale point-to-multipole transform */
  if (reg->pm != NULL) {
    fxt_matld_scalecol(reg->pm, ss);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* scale point-to-point and point-to-local transforms */
  for (cell = reg->exp_list; cell != NULL; cell = cell->exp_next)
    if (cell->type == FMM_CTLOC || cell->type == FMM_CTDIR) {
      fxt_matld_scalecol(cell->m, ss);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

  /* deallocate vector */
  fxt_vecld_del(ss);

  /* recursion on children */
  if (reg->ch0 != NULL) {
    scalecol_reg(reg->ch0, s);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
    scalecol_reg(reg->ch1, s);
  }
}


/*** scale columns ***/
void fxt_fmmld_scalecol(fxt_fmmld *fmm, fxt_vecld *s) {
  int i;

  /* check null pointers */
  if (fmm == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_scalecol: null pointer\n");
    return;
  }

  /* check size */
  if (fmm->ns != s->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_scalecol: size mismatch\n");
    return;
  }

  /* scale original matrix */
  fxt_matld_scalecol(fmm->mat, s);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* scale FMM */
  for (i=0; i< fmm->n_roots; i++) {
    scalecol_reg(&fmm->roots[i], s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}


/*** scale target points ***/
static void scalerow_reg(fmmld_reg *reg, fxt_vecld *s) {
  fmmld_cell *cell;
  fxt_vecld *ss;

  if (reg->nt == 0)
    return;

  /* get the corresponding part of scaling vector */
  ss = fxt_vecld_clone(s, reg->it, reg->it + reg->nt - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* scale local to point transform matrix */
  if (reg->lp != NULL) {
    fxt_matld_scalerow(reg->lp, ss);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* point-to-point and multipole-to-point transforms */
  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next)
    if (cell->type == FMM_CTDIR || cell->type == FMM_CTMUL) {
      fxt_matld_scalerow(cell->m, ss);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

  /* deallocate vector */
  fxt_vecld_del(ss);

  /* recursion to children */
  if (reg->ch0 != NULL) {
    scalerow_reg(reg->ch0, s);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
    scalerow_reg(reg->ch1, s);
  }
}


/*** scale rows of FMM ***/
void fxt_fmmld_scalerow(fxt_fmmld *fmm, fxt_vecld *s) {
  int i;

  /* check null pointers */
  if (fmm == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_scalerow: null pointer\n");
    return;
  }

  /* check size */
  if (fmm->nt != s->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_scalerow: size mismatch\n");
    return;
  }

  /* scale original matrix */
  fxt_matld_scalerow(fmm->mat, s);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* scale FMM */
  for (i=0; i< fmm->n_roots; i++) {
    scalerow_reg(&fmm->roots[i], s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}
