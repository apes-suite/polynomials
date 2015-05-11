#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  expansion of FMM
*********************************************************************/

/*** generate local expansions ***/
static int exp_loc(fmmld_reg *reg, fxt_vecld *u,
		   fmmld_reg *preg, int pcnt) {
  fmmld_reg *lreg;		/* region of local expansion */
  int lcnt;			/* accumulation counter */

  if (reg->type & FMM_RTLOC) {
    
    /* local expansion of this region */
    lreg = reg;
    lcnt = 0;
  } else {

    /* local expansion of parent region */
    lreg = preg;
    lcnt = pcnt;
  }

  if ((reg->type & FMM_RTEVL) && lreg != NULL) {

    /* point to local expansion */
    if ((lcnt++) == 0)
      fxt_matld_settrmatprtvec(lreg->cl, reg->lp,
			       u, reg->it, reg->it + reg->nt - 1);
    else
      fxt_matld_addtrmatprtvec(lreg->cl, reg->lp,
			       u, reg->it, reg->it + reg->nt - 1);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    lreg = NULL;
  } 

  if (reg->ch0 != NULL) {

    /* recursion to children */
    lcnt = exp_loc(reg->ch0, u, lreg, lcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    lcnt = exp_loc(reg->ch1, u, lreg, lcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0;
  }

  /* local to local transform */
  if ((reg->type & FMM_RTLOC) && preg != NULL) {

    /* set parent's local expansion */
    if ((pcnt++) == 0)
      fxt_matld_settrmatvec(preg->cl, reg->ll, reg->cl);
    else
      fxt_matld_addtrmatvec(preg->cl, reg->ll, reg->cl);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* parent's accumulation counter */
    lcnt = pcnt;
  }

  return lcnt;
}
  

/*** create multipole expansions and evaluate them ***/
static void exp_mul(fmmld_reg *reg, fxt_vecld *u, fxt_vecld *v,
		    fmmld_reg *preg, int pcnt) {
  fmmld_cell *cell;
  int mcnt = 0;			/* accumulation counter */

  /* create multipole expansion */
  if (reg->type & FMM_RTMUL) {

    if (preg != NULL) {
      /* this must be the first accumulation */
      mcnt ++;

      /* multipole-to-multipole transform */
      fxt_matld_settrmatvec(reg->cm, reg->mm, preg->cm);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    /* multipole expansion of this region */
    preg = reg;
  }

  /* transform and/or evaluate from other regions */
  for (cell = reg->exp_list; cell != NULL; cell = cell->exp_next) {
    fmmld_reg *creg = cell->exp_reg;

    if (cell->type == FMM_CTTRN) {

      /* local to multipole transform */
      if ((mcnt++) == 0)
	fxt_matld_settrmatvec(reg->cm, cell->m, creg->cl);
      else
	fxt_matld_addtrmatvec(reg->cm, cell->m, creg->cl);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTMUL) {

      /* point to multipole expansion */
      if ((mcnt++) == 0)
	fxt_matld_settrmatprtvec(reg->cm, cell->m, u,
				 creg->it, creg->it + creg->nt - 1);
      else
	fxt_matld_addtrmatprtvec(reg->cm, cell->m, u,
				 creg->it, creg->it + creg->nt - 1);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTDIR) {

      /* point to point evaluation */
      if ((pcnt++) == 0)
	fxt_matld_prtsettrmatprtvec(v, reg->is, reg->is + reg->ns - 1,
				    cell->m, u, creg->it,
				    creg->it + creg->nt - 1);
      else
	fxt_matld_prtaddtrmatprtvec(v, reg->is, reg->is + reg->ns - 1,
				    cell->m, u, creg->it,
				    creg->it + creg->nt - 1);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTLOC) {

      /* direct evaluation of local expansion */
      if ((pcnt++) == 0)
	fxt_matld_prtsettrmatvec(v, reg->is, reg->is + reg->ns - 1,
				 cell->m, creg->cl);
      else
	fxt_matld_prtaddtrmatvec(v, reg->is, reg->is + reg->ns - 1,
				 cell->m, creg->cl);

      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }
  }

  /* evaluate multipole expansion */
  if ((reg->type & FMM_RTEXP) && preg != NULL) {

    /* evaluate multipole expansion */
    if ((pcnt++) == 0)
      fxt_matld_prtsettrmatvec(v, reg->is, reg->is + reg->ns - 1,
			       reg->pm, preg->cm);
    else
      fxt_matld_prtaddtrmatvec(v, reg->is, reg->is + reg->ns - 1,
			       reg->pm, preg->cm);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    preg = NULL;
  }

  if (reg->ch0 != NULL) {

    /* recursion to children */
    exp_mul(reg->ch0, u, v, preg, pcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    exp_mul(reg->ch1, u, v, preg, pcnt);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}


/*** expansion (transposed evaluation) of FMM ***/
void fxt_fmmld_exp(fxt_vecld *v, fxt_fmmld *fmm, fxt_vecld *u) {
  int i;

  /* check null pointers */
  if (v == NULL || fmm == NULL || u == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_exp: null pointer\n");
    return;
  }

  /* check size */
  if (fmm->ns != v->n || fmm->nt != u->n) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_exp: size mismatch\n");
    return;
  }

  /* generate local expansions */
  for (i=0; i< fmm->n_roots; i++) {
    exp_loc(&fmm->roots[i], u, NULL, 0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* generate multipole expansions and evaluate them */
  for (i=0; i< fmm->n_roots; i++) {
    exp_mul(&fmm->roots[i], u, v, NULL, 0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}
