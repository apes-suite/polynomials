#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generate multipole expansion matrices for FMM
*********************************************************************/


/*** computes far-field for the multipole expansion ***/
static void get_farfield(fmmld_reg *reg, long x[2]) {
  fmmld_cell *cell;

  /* get the father's far field */
  if (reg->preg != NULL)
    get_farfield(reg->preg, x);
  else {
    x[0] = 0;  x[1] = reg->pfmm->nt;
  }

  /* check transform matrices */
  for (cell = reg->exp_list; cell != NULL; cell = cell->exp_next)
    if (cell->type == FMM_CTTRN || cell->type == FMM_CTMUL) {
      fmmld_reg *creg = cell->exp_reg;

      if (creg->it < reg->it) {	/* left side */
	if (x[0] < creg->it + creg->nt)
	  x[0] = creg->it + creg->nt;
      } else {			/* right side */
	if (x[1] > creg->it)
	  x[1] = creg->it;
      }
    }
}


/*** get size of the multipole expansion ***/
static long get_mulsize(fmmld_reg *reg) {

  /* no point, no expansion */
  if (reg->ns == 0)
    return 0;

  /* having multipole expansion */
  if (reg->type & FMM_RTMUL)
    return reg->km;

  /* no expansion, but an expanding region */
  if (reg->type & FMM_RTEXP)
    return reg->ns;

  /* resort to children */
  return get_mulsize(reg->ch0) + get_mulsize(reg->ch1);
}


/*** get the multipole evaluation matrix ***/
static long get_mulmat(fmmld_reg *reg, long n, fxt_matld *a) {

  /* no point, no expansion */
  if (reg->ns == 0)
    return n;

  /* having multipole expansion */
  if (reg->type & FMM_RTMUL) {
    fxt_vecld *s;  fxt_matld *p;

    /* void expansion */
    if (reg->km == 0)
      return n;

    /* allocate matrix */
    p = fxt_matld_clone(reg->mp, 0, reg->pfmm->nt - 1,
			0, reg->km - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* get scaling vector */
    s = fxt_fmmld_get_scalemul(reg);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* scaling */
    fxt_matld_scalecol(p, s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* put matrix to A */
    fxt_matld_prtset(a, 0, reg->pfmm->nt - 1, n, n + reg->km - 1, p);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* deallocate temporals */
    fxt_vecld_del(s);
    fxt_matld_del(p);

    /* return next point */
    return n + reg->km;
  }

  /* no expansion, but an expanding region */
  if (reg->type & FMM_RTEXP) {
    /* copy original matrix */
    fxt_matld_prtsetprt(a, 0, reg->pfmm->nt - 1, n, n + reg->ns - 1,
			reg->pfmm->mat, 0, reg->pfmm->nt - 1,
			reg->is, reg->is + reg->ns - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* return next point */
    return n + reg->ns;
  }

  /* resort to children */
  n = get_mulmat(reg->ch0, n, a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  n = get_mulmat(reg->ch1, n, a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  return n;
}


/*** put the multipole transform matrix ***/
static long put_mulmat(fmmld_reg *reg, long n, fxt_matld *a) {

  /* no point, no expansion */
  if (reg->ns == 0)
    return n;

  /* having multipole expansion */
  if (reg->type & FMM_RTMUL) {

    /* first, resize matrix */
    fxt_matld_shrink(reg->mm, a->nrow, reg->km);

    if (a->nrow > 0 && reg->km > 0) {
      /* get the scaling vector */
      fxt_vecld *s = fxt_fmmld_get_scalemul(reg);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0;

      /* inversion */
      fxt_vecld_einv(s);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0;

      /* copy matrix */
      fxt_matld_setprt(reg->mm,
		       a, 0, a->nrow - 1, n, n + reg->km - 1);

      /* re-scale matrix */
      fxt_matld_scalecol(reg->mm, s);

      /* deallocate vector */
      fxt_vecld_del(s);
    }

    /* return next point */
    return n + reg->km;
  }

  /* no expansion, but an expanding region */
  if (reg->type & FMM_RTEXP) {

    /* resize matrix */
    fxt_matld_shrink(reg->pm, a->nrow, reg->ns);

    /* copy matrix */
    if (a->nrow > 0 && reg->ns > 0)
      fxt_matld_setprt(reg->pm,
		       a, 0, a->nrow - 1, n, n + reg->ns - 1);

    /* return next point */
    return n + reg->ns;
  }

  /* resort to children */
  n = put_mulmat(reg->ch0, n, a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  n = put_mulmat(reg->ch1, n, a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* return next point */
  return n;
}


/*** make multipole expansion matrices ***/
void fxt_fmmld_make_mulexpmat(fmmld_reg *reg) {

  /* check null pointers */
  if (reg == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_make_mulexpmat: null pointer\n");
    return;
  }

  /* recursive call to children */
  if (reg->ch0 != NULL) {
    fxt_fmmld_make_mulexpmat(reg->ch0);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    fxt_fmmld_make_mulexpmat(reg->ch1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  if (reg->type & FMM_RTMUL) {
    /* create multipole expansion */

    if (reg->type & FMM_RTEXP) {
      /* point to multipole expansion creation */

      long x[2];  fxt_matld *a;
      int errc;

      /* get the far-field points */
      get_farfield(reg, x);

      /* create matrix */
      a = fxt_matld_new(reg->pfmm->nt, reg->ns);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* set first far-field matrices */
      if (0 < x[0]) {
	fxt_matld_prtsetprt(a, 0, x[0] - 1, 0, reg->ns - 1,
			    reg->pfmm->mat, 0, x[0] - 1,
			    reg->is, reg->is + reg->ns - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

      /* clear near-field region */
      if (x[0] < x[1]) {
	fxt_matld_zero(a, x[0], x[1] - 1, 0, reg->ns - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

      /* set second far-field matrices */
      if (x[1] < reg->pfmm->nt) {
	fxt_matld_prtsetprt(a, x[1], reg->pfmm->nt - 1, 0, reg->ns -1,
			    reg->pfmm->mat, x[1], reg->pfmm->nt - 1,
			    reg->is, reg->is + reg->ns - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }

      /* get expansion and evaluation matrices */
      errc = fxt_matld_lowrank(a, reg->pfmm->eps,
			       reg->mp, reg->pm);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      /* mark insufficient multipole expansion order */
      if (errc == 1)
	reg->errc |= 1;

      reg->km = reg->mp->ncol;

      /* deallocate matrix */
      fxt_matld_del(a);
    } else {			/* reg->type & FMM_RTEXP == 0 */
      /* multipole expansion by transform */

      /* get size of multipole expansions */
      long m = get_mulsize(reg->ch0) + get_mulsize(reg->ch1);

      if (m > 0) {		/* the normal case */
	fxt_matld *a, *z;
	long n, x[2];
	int errc;

	/* allocate matrix */
	a = fxt_matld_new(reg->pfmm->nt, m);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get the multipole evaluation matrix of ch0 */
	n = get_mulmat(reg->ch0, 0, a);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get the multipole evaluation matrix of ch1 */
	n = get_mulmat(reg->ch1, n, a);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get the near-field region */
	get_farfield(reg, x);

	/* clear near-field region */
	if (x[0] < x[1]) {
	  fxt_matld_zero(a, x[0], x[1] - 1, 0, m - 1);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;
	}

	/* allocate the transform matrix */
	z = fxt_matld_new(reg->km, m);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* get transform and evaluation matrices */
	errc = fxt_matld_lowrank(a, reg->pfmm->eps, reg->mp, z);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* mark insufficient multipole order */
	if (errc == 1)
	  reg->errc |= 1;

	reg->km = reg->mp->ncol;

	/* put the multipole transform matrix of ch0 */
	n = put_mulmat(reg->ch0, 0, z);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* put the multipole transform matrix of ch1 */
	n = put_mulmat(reg->ch1, n, z);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;

	/* deallocate temporal matrices */
	fxt_matld_del(z);
	fxt_matld_del(a);
      } else {			/* m == 0: ignorable submatrix */

	fxt_matld *z = fxt_matld_new(0, 0);

	/* put void for multipole transform */
	put_mulmat(reg->ch0, 0, z);
	put_mulmat(reg->ch1, 0, z);

	fxt_matld_del(z);

	/* zero-sized matrix for multipole evaluation */
	fxt_matld_shrink(reg->mp, reg->pfmm->nt, 0);
	reg->km = 0;

	fxt_error_raise();
      }
    }
  } /* if (reg->type & FMM_RTMUL) */
}
