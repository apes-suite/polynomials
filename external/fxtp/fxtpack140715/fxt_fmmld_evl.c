#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  evaluation of FMM
*********************************************************************/

/*** compute multipole expansions ***/
static int evl_mul(fmmld_reg *reg, fxt_vecld *u,
		   fmmld_reg *preg, int pcnt) {
  fmmld_reg *mreg;		/* region for multipole expansion */
  int mcnt;			/* accumulation count for mul-exp */

  if (reg->type & FMM_RTMUL) {

    /* expansion of this region */
    mreg = reg;
    mcnt = 0;
  } else {

    /* expansion of a parent region */
    mreg = preg;
    mcnt = pcnt;
  }

  if ((reg->type & FMM_RTEXP) && mreg != NULL) {

    /* generate multipole expansion */
    if ((mcnt++) == 0)
      fxt_matld_setmatprtvec(mreg->cm, reg->pm,
			     u, reg->is, reg->is + reg->ns - 1);
    else
      fxt_matld_addmatprtvec(mreg->cm, reg->pm,
			     u, reg->is, reg->is + reg->ns - 1);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* expansion is done */
    mreg = NULL;
  }

  if (reg->ch0 != NULL) {

    /* recursion to children */
    mcnt = evl_mul(reg->ch0, u, mreg, mcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    mcnt = evl_mul(reg->ch1, u, mreg, mcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;
  }

  /* transform to parent's multipole expansion */
  if ((reg->type & FMM_RTMUL) && preg != NULL) {

    /* transform to parent's multipole expansion */
    if ((pcnt++) == 0)
      fxt_matld_setmatvec(preg->cm, reg->mm, reg->cm);
    else
      fxt_matld_addmatvec(preg->cm, reg->mm, reg->cm);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    mcnt = pcnt;
  }

  /* parent's accumulation counter */
  return mcnt;
}


/*** compute local expansion and evaluation ***/
static void evl_loc(fmmld_reg *reg, fxt_vecld *u, fxt_vecld *v,
		    fmmld_reg *preg, int pcnt) {
  fmmld_cell *cell;  int lcnt;

  /* local expansion accumulation counter */
  lcnt = 0;

  if (reg->type & FMM_RTLOC) {
    if (preg != NULL) {
      /* this must be the first accumulation */
      lcnt++;

      /* transorm from parent's local expansion */
      fxt_matld_setmatvec(reg->cl, reg->ll, preg->cl);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    /* having local expansion */
    preg = reg;
  }

  /* transform and/or evaluate from other regions */
  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;

    if (cell->type == FMM_CTTRN) {

      /* multipole-to-local transform */
      if ((lcnt++) == 0)
	fxt_matld_setmatvec(reg->cl, cell->m, creg->cm);
      else
	fxt_matld_addmatvec(reg->cl, cell->m, creg->cm);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTLOC) {

      /* point to local expansion */
      if ((lcnt++) == 0)
	fxt_matld_setmatprtvec(reg->cl, cell->m,
			       u, creg->is, creg->is + creg->ns - 1);
      else
	fxt_matld_addmatprtvec(reg->cl, cell->m,
			       u, creg->is, creg->is + creg->ns - 1);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTDIR) {

      /* point to point evaluation */
      if ((pcnt++) == 0)
	fxt_matld_prtsetmatprtvec(v, reg->it, reg->it + reg->nt - 1,
				  cell->m, u,
				  creg->is, creg->is + creg->ns - 1);
      else
	fxt_matld_prtaddmatprtvec(v, reg->it, reg->it + reg->nt - 1,
				  cell->m, u,
				  creg->is, creg->is + creg->ns - 1);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTMUL) {

      /* direct evaluation of multipole expansion */
      if ((pcnt++) == 0)
	fxt_matld_prtsetmatvec(v, reg->it, reg->it + reg->nt - 1,
			       cell->m, creg->cm);
      else
	fxt_matld_prtaddmatvec(v, reg->it, reg->it + reg->nt - 1,
			       cell->m, creg->cm);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }
  }

  if ((reg->type & FMM_RTEVL) && preg != NULL) {

    /* evaluate the accumulated local expansion */
    if ((pcnt++) == 0)
      fxt_matld_prtsetmatvec(v, reg->it, reg->it + reg->nt - 1,
			     reg->lp, preg->cl);
    else
      fxt_matld_prtaddmatvec(v, reg->it, reg->it + reg->nt - 1,
			     reg->lp, preg->cl);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    preg = NULL;
  }

  if (reg->ch0 != NULL) {

    /* recursion to children */
    evl_loc(reg->ch0, u, v, preg, pcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    evl_loc(reg->ch1, u, v, preg, pcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}


/*** evaluate FMM ***/
void fxt_fmmld_evl(fxt_vecld *v, fxt_fmmld *fmm, fxt_vecld *u) {
  int i;

  /* check null pointers */
  if (v == NULL || fmm == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_evl: null pointer\n");
    return;
  }

  /* check size */
  if (fmm->ns != u->n || fmm->nt != v->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_evl: size mismatch\n");
    return;
  }

  /* generate multipole expansion */
  for (i=0; i< fmm->n_roots; i++) {
    evl_mul(&fmm->roots[i], u, NULL, 0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* generate local expansion and evaluate them */
  for (i=0; i< fmm->n_roots; i++) {
    evl_loc(&fmm->roots[i], u, v, NULL, 0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}
