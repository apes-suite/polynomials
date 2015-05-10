#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  preprocessing driver routines for FMM
*********************************************************************/

#define MAX(x, y) ((x) < (y) ? (x) : (y))

#define UP_RATE   1.78
#define DOWN_RATE 0.56
#define EPSILON   1e-16

static void opt_eps(fxt_fmmld*);

/*** preprocessing with optimized parameter ***/
void fxt_fmmld_preproc(fxt_fmmld *fmm) {
  double e, e2;			/* the error */

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_preproc: null pointer\n");
    return;
  }

  /* initial precision parameter */
  fmm->eps = fmm->prec * UP_RATE;

  /* initial preprocessing with that eps */
  opt_eps(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* compute the error norm of the resulting FMM */
  e = fxt_fmmld_error(fmm, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* for safety, recheck error */
  e2 = fxt_fmmld_error(fmm, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (e < e2) e = e2;

#if 0
  /* increase eps until error >= required precision */
  while (e < fmm->prec && fmm->eps < fmm->cnorm
	 && fmm->eps < fmm->prec * 100.0) {

    /* modify eps according to the result */
    if (e < EPSILON * fmm->cnorm)
      fmm->eps *= UP_RATE;
    else
      fmm->eps *= fmm->prec / e * UP_RATE;

    /* clear previous data structure */
    fxt_fmmld_unpreproc(fmm);

    /* remake with the new eps */
    opt_eps(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* compute the error */
    e = fxt_fmmld_error(fmm, POWERPREC);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }
#endif

  /* decrease eps until error <= required precision */
  while (e * (1.0 + POWERPREC) > fmm->prec
	 && fmm->eps > fmm->prec * 0.01) {
    fmm->eps *= DOWN_RATE;

    /* clear previous data structure */
    fxt_fmmld_unpreproc(fmm);

    /* remake with the new eps */
    opt_eps(fmm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* compute the error */
    e = fxt_fmmld_error(fmm, POWERPREC);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* recompute error for safety */
    e2 = fxt_fmmld_error(fmm, POWERPREC);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    if (e < e2) e = e2;
  }

  /* once unmake matrix */
  fxt_fmmld_unmake_matrix(fmm);

  /* remake matrix with true memory requirements */
  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (e > fmm->prec)
    fmm->warn = 1;
}


/*** undo preprocessing ***/
void fxt_fmmld_unpreproc(fxt_fmmld *fmm) {

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_unpreproc: null pointer\n");
    return;
  }

  fxt_fmmld_unmake_matrix(fmm);
  fxt_fmmld_unmake_network(fmm);
  fxt_fmmld_unmake_regions(fmm);
}


/*** make estimate ***/
static long make_estimate(fxt_fmmld *fmm) {
  fxt_fmmld_make_regions(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  fxt_fmmld_make_network(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  return fxt_fmmld_estimate_mflop(fmm);
}


/*** undo estimate ***/
static void unmake_estimate(fxt_fmmld *fmm) {
  fxt_fmmld_unmake_network(fmm);
  fxt_fmmld_unmake_regions(fmm);
}


static long best_flop;
static int best_apptype;
static int best_nroots;
static int best_depth;


/*** save the current parameters ***/
static void save_best(fxt_fmmld *fmm, long flop) {
  best_flop = flop;
  best_apptype = fmm->approx_type;
  best_nroots = fmm->n_roots;
  best_depth = fmm->max_depth;
}


/*** read the saved parameters ***/
static void read_best(fxt_fmmld *fmm) {
  fmm->approx_type = best_apptype;
  fmm->n_roots = best_nroots;
  fmm->max_depth = best_depth;
}


/*** search for best parameter ***/
static void search_param(fxt_fmmld *fmm, int gen_flag) {
  int apptype, nmax, n_roots, dmax, depth;
  long np, flop;  double looseness;
  static long flop_bound = 0;

  best_flop = -1;

  if (fmm->kest < 2)
    looseness = 3.0;
  else
    looseness = (fmm->kest + 1.0) / (fmm->kest - 1.0);

  looseness *= looseness;

  /* try all approximation types */
  for (apptype = 1; apptype <= 3; apptype ++) {
    fmm->approx_type = apptype;

    /* compute maximum for n_roots */
    nmax = (fmm->ns + fmm->nt) / fmm->min_np / 2;
    if (nmax < 3)
      nmax = 3;

    /* try all numbers of roots */
    for (n_roots = 3; n_roots <= nmax; n_roots ++) {
      fmm->n_roots = n_roots;

      /* compute maximum tree depth */
      np = MAX(fmm->ns, fmm->nt) / fmm->n_roots;
      for (dmax = 0; np >= fmm->min_np; dmax++)
	np /= 2;

      /* try all depth limits */
      for (depth = 0; depth <= dmax; depth ++) {
	fmm->max_depth = depth;

	/* try making FMM */
	if (gen_flag == 0) {
	  /* make rough estimate */
	  flop = make_estimate(fmm);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* undo estimate */
	  unmake_estimate(fmm);

	  /* save if it is the best */
	  if (best_flop < 0 || flop < best_flop)
	    save_best(fmm, flop);
	} else {
	  /* make rough estimate */
	  flop = make_estimate(fmm);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return;

	  /* make if promisful */
	  if (flop <= flop_bound * looseness) {

	    /* make matrix */
	    fxt_fmmld_make_matrix(fmm);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* evaluate flops */
	    flop = fxt_fmmld_evaluate_mflop(fmm);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return;

	    /* save if it is the best */
	    if (best_flop < 0 || flop < best_flop)
	      save_best(fmm, flop);

	    /* undo preprocessing */
	    fxt_fmmld_unmake_matrix(fmm);
	  }

	  /* undo estimate */
	  unmake_estimate(fmm);
	}
      }
    }
  }

  flop_bound = best_flop;
}
    

/*** optimum preprocessing with the given eps ***/
static void opt_eps(fxt_fmmld *fmm) {
  int d;  long n;

  /* estimate kmax */
  fxt_fmmld_estimate_kmax(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* set default parameters */
  fmm->approx_type = 3;
  fmm->n_roots = 3;

  n = MAX(fmm->ns, fmm->nt) / 3;
  for (d = 0; n >= fmm->min_np; d++, n /= 2);
  fmm->max_depth = d;

  fmm->min_np = 2;

  /* test run */
  fxt_fmmld_make_regions(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_network(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* evaluate the default make */
  fxt_fmmld_evaluate_kmax(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  if (fmm->kest < 0)
    return;

  fmm->kmax = (fmm->kmax + 1) * 2;
  if (fmm->kmax < 6)
    fmm->kmax = 6;

  /* undo make */
  fxt_fmmld_unmake_matrix(fmm);
  fxt_fmmld_unmake_network(fmm);
  fxt_fmmld_unmake_regions(fmm);

  /* make rough estimate */
  search_param(fmm, 0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* try preprocessing and compare */
  search_param(fmm, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* make the best one */
  read_best(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* make the regions */
  fxt_fmmld_make_regions(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* make the network */
  fxt_fmmld_make_network(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* make the matrix */
  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  return;
}
