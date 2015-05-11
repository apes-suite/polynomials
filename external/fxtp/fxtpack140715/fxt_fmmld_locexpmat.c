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

  generate local expansion matrices for FMM
*********************************************************************/


/*** sparsify matrix ***/
static void sparsify(fxt_matld *a, double eps) {
  double th;

  th = fxt_matld_dropthreshould(a, eps);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_matld_drop(a, th);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}


/*** make local expansions ***/
void fxt_fmmld_make_locexpmat(fmmld_reg *reg) {
  fmmld_reg *preg;
  fmmld_cell *cell;

  /* check null pointers */
  if (reg == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_make_locexpmat: null pointer\n");
    return;
  }

  /* search for the lowest predecessor with local expansion */
  for (preg = reg->preg; preg != NULL; preg = preg->preg)
    if (preg->type & FMM_RTEVL) {
      preg = NULL;
      break;
    } else if (preg->type & FMM_RTLOC)
      break;

  if (reg->type & FMM_RTLOC) {
    /* generate local expansion */

    long m;			/* matrix size */

    m = (preg != NULL ? preg->kl : 0);

    for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next)
      if (cell->type == FMM_CTTRN)
	m += cell->evl_reg->km;
      else if (cell->type == FMM_CTLOC)
	m += cell->evl_reg->ns;

    if (m > 0) {		/* normal case */
      fxt_matld *a, *z;  long n;  int errc;

      /* allocate matrix */
      a = fxt_matld_new(reg->nt, m);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* collect matrix from preg */
      if (preg != NULL && preg->kl > 0) {
	fxt_matld *p;
	fxt_vecld *s;

	/* cut from parent's evaluation matrix */
	p = fxt_matld_clone(preg->lp, reg->it - preg->it,
			    reg->it - preg->it + reg->nt - 1,
			    0, preg->kl - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get the scaling vector */
	s = fxt_fmmld_get_scaleloc(preg);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* scale the matrix */
	fxt_matld_scalecol(p, s);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* store the results in A */
	fxt_matld_prtset(a, 0, reg->nt - 1, 0, preg->kl - 1, p);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
	
	/* deallocate temporals */
	fxt_vecld_del(s);
	fxt_matld_del(p);

	/* the next column */
	n = preg->kl;
      } else
	n = 0;
      
      /* collect matrix from transform lists */
      for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next){
	fmmld_reg *creg = cell->evl_reg;

	if (cell->type == FMM_CTTRN && creg->km > 0) {
	  fxt_matld *p;
	  fxt_vecld *s;

	  /* cut from creg's evaluation matrix */
	  p = fxt_matld_clone(creg->mp,
			      reg->it, reg->it + reg->nt - 1,
			      0, creg->km - 1);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;
	  
	  /* get the scaling vector */
	  s = fxt_fmmld_get_scalemul(creg);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* scale the matrix */
	  fxt_matld_scalecol(p, s);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* store the results in A */
	  fxt_matld_prtset(a, 0, reg->nt - 1, n, n + creg->km - 1, p);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;
	  
	  /* deallocate temporals */
	  fxt_vecld_del(s);
	  fxt_matld_del(p);

	  /* next column */
	  n += creg->km;
	} else if (cell->type == FMM_CTLOC) {
	  /* just copy from original matrix */
	  fxt_matld_prtsetprt(a, 0, reg->nt - 1, n, n + creg->ns - 1,
			      reg->pfmm->mat,
			      reg->it, reg->it + reg->nt - 1,
			      creg->is, creg->is + creg->ns - 1);

	  n += creg->ns;
	}
      }

      /* allocate temporal creation matrix */
      z = fxt_matld_new(reg->kl, n);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* make creation and evaluation matrices */
      errc = fxt_matld_lowrank(a, reg->pfmm->eps, reg->lp, z);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* mark insufficient local expansion order */
      if (errc == 1)
	reg->errc |= 2;

      reg->kl = reg->lp->ncol;

      /* put local-local transform matrix */
      if (preg != NULL) {
	/* resize matrix */
	fxt_matld_shrink(reg->ll, reg->kl, preg->kl);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	if (preg->kl > 0 && reg->kl > 0) {
	  fxt_vecld *s;

	  /* set local-local transform matrix  */
	  fxt_matld_setprt(reg->ll,
			   z, 0, reg->kl - 1, 0, preg->kl - 1);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* get the scaling vector */
	  s = fxt_fmmld_get_scaleloc(preg);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* inversion for re-scaling */
	  fxt_vecld_einv(s);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* descale the matrix */
	  fxt_matld_scalecol(reg->ll, s);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* deallocate temporal */
	  fxt_vecld_del(s);
	}

	/* next column */
	n = preg->kl;
      } else
	n = 0;

      /* put matrix to transform lists */
      for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next){
	fmmld_reg *creg = cell->evl_reg;

	if (cell->type == FMM_CTTRN) {

	  /* resize matrix */
	  fxt_matld_shrink(cell->m, reg->kl, creg->km);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  if (reg->kl > 0 && creg->km > 0) {
	    fxt_vecld *s;

	    /* get the transform matrix */
	    fxt_matld_setprt(cell->m,
			     z, 0, reg->kl - 1, n, n + creg->km - 1);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* get the scaling vector */
	    s = fxt_fmmld_get_scalemul(creg);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* inversion for descaling */
	    fxt_vecld_einv(s);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* descale the matrix */
	    fxt_matld_scalecol(cell->m, s);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* deallocate temporal */
	    fxt_vecld_del(s);
	  }

	  /* next column */
	  n += creg->km;
	} else if (cell->type == FMM_CTLOC) {

	  /* resize matrix */
	  fxt_matld_shrink(cell->m, reg->kl, creg->ns);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  if (reg->kl > 0) {
	    /* set generation matrix */
	    fxt_matld_setprt(cell->m,
			     z, 0, reg->kl - 1, n, n + creg->ns - 1);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;
	  }

	  /* next column */
	  n += creg->ns;
	}
      }

      /* deallocate temporals */
      fxt_matld_del(z);
      fxt_matld_del(a);
    } else {			/* m == 0: ignorable submatrix */

      /* zero-sized expansion */
      reg->kl = 0;

      /* resize matrix */
      fxt_matld_shrink(reg->lp, reg->nt, reg->kl);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      if (preg != NULL) {
	/* resize matrix */
	fxt_matld_shrink(reg->ll, reg->kl, preg->kl);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

      for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next){
	fmmld_reg *creg = cell->evl_reg;

	if (cell->type == FMM_CTTRN) {

	  /* resize matrix */
	  fxt_matld_shrink(cell->m, reg->kl, creg->km);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	} else if (cell->type == FMM_CTLOC) {

	  /* resize matrix */
	  fxt_matld_shrink(cell->m, reg->kl, creg->ns);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;
	}
      }
    }

    /* I have local expansion */
    preg = reg;
  } /* if (reg->type & FMM_RTLOC */

  /* create evaluation matrix */
  if (preg != NULL && (reg->type & FMM_RTEVL)) {

    /* resize matrix */
    fxt_matld_shrink(reg->lp, reg->nt, preg->kl);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    if (reg->nt > 0 && preg->kl > 0) {
      fxt_vecld *s;

      if (reg != preg) {
	/* get the matrix */
	fxt_matld_setprt(reg->lp, preg->lp, reg->it - preg->it,
			 reg->it - preg->it + reg->nt - 1,
			 0, preg->kl - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

      /* get the scaling vector */
      s = fxt_fmmld_get_scaleloc(preg);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* scale the matrix */
      fxt_matld_scalecol(reg->lp, s);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* sparsify the matrix */
      sparsify(reg->lp, reg->pfmm->eps);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* inversion for descaling */
      fxt_vecld_einv(s);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* descale the matrix */
      fxt_matld_scalecol(reg->lp, s);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
      
      /* deallocate temporal */
      fxt_vecld_del(s);
    }
  }

  /* generate evaluation matrices */
  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;

    if (cell->type == FMM_CTDIR) {

      /* copy from the original matrix */
      fxt_matld_setprt(cell->m, reg->pfmm->mat,
		       reg->it, reg->it + reg->nt - 1,
		       creg->is, creg->is + creg->ns - 1);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* sparsify if possible */
      sparsify(cell->m, reg->pfmm->eps);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

    } else if (cell->type == FMM_CTMUL) {

      /* resize matrix */
      fxt_matld_shrink(cell->m, reg->nt, creg->km);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      if (reg->nt > 0 && creg->km > 0) {
	fxt_vecld *s;

	/* get the evaluation matrix */
	fxt_matld_setprt(cell->m,
			 creg->mp, reg->it, reg->it + reg->nt - 1,
			 0, creg->km - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get the scaling vector */
	s = fxt_fmmld_get_scalemul(creg);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* scale matrix */
	fxt_matld_scalecol(cell->m, s);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* sparsify matrix */
	sparsify(cell->m, reg->pfmm->eps);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* inversion for descale */
	fxt_vecld_einv(s);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* descaling matrix */
	fxt_matld_scalecol(cell->m, s);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* deallocate temporal */
	fxt_vecld_del(s);
      }
    }
  }

  /* recursion to children */
  if (reg->ch0 != NULL) {
    fxt_fmmld_make_locexpmat(reg->ch0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_fmmld_make_locexpmat(reg->ch1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
}

