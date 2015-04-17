#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  make FXT structure
*********************************************************************/

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MINSIZE 10

#if SHOWFXTMAKE
static void print_depth(fxtld *fxt) {
  if (fxt->par == NULL)
    printf(".");
  else {
    print_depth(fxt->par);

    if (fxt->par->ch0 == fxt)
      printf("L");
    else
      printf("H");
  }
}
#endif

static fxt_matll* get_scaled_bmat_ll(fxtld *fxt) {
  fxt_matll *mat;
  fxt_vecll *a1;
  fxt_vecld *a0;

  fxt_matll_permrow(fxt->mat, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  mat = fxt_matll_clone(fxt->mat, 0, fxt->mat->ncol - 1,
			0, fxt->mat->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matll_ipermrow(fxt->mat, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  a0 = fxtld_get_ascale(fxt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  a1 = fxt_vecld_to_vecll(a0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matll_scalerow(mat, a1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecll_del(a1);
  fxt_vecld_del(a0);

  return mat;
}

static fxt_matld* get_scaled_bmat_ld(fxtld *fxt) {
  fxt_matld *mat, *mat0;
  fxt_vecld *a;

  mat = fxt_matld_new(fxt->mat->ncol, fxt->mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  mat0 = fxt_matll_to_matld(fxt->mat);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_permrow(mat0, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_setprt(mat, mat0, 0, mat->ncol - 1, 0, mat->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_del(mat0);

  a = fxtld_get_ascale(fxt);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_matld_scalerow(mat, a);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  fxt_vecld_del(a);

  return mat;
}

long fxtld_make(fxtld *fxt, fxt_matll *mat, fxt_vecll *p,
		double eps, double th, long flop0, long ui0) {
  long flop_c, flop_r, flop_d, flop_i, flop_id, flop_iv;
  
  /* check null pointers */
  if (fxt == NULL || mat == NULL || p == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_make: null pointer\n");
    return 0;
  }

  /* range of input vector */
  fxt->ui0 = ui0;
  fxt->ui1 = ui0 + mat->ncol - 1;
  fxt->mat = mat;
  fxt->p = p;

  /* initialize fields */
  fxt->n_drop = 0;
  fxt->idx = NULL;
  fxt->fmm = NULL;
  fxt->imat = NULL;
  fxt->xs = fxt->xt = NULL;

  fxt->de = NULL;
  fxt->ch0 = fxt->ch1 = NULL;
  fxt->bmat = NULL;
  fxt->pp = NULL;
  fxt->ws = fxt->wt = NULL;

  /* check direct evaluation */
  flop_r = flop_d = fxt_matll_dropflop(mat, th);
  flop_c = (flop0 < flop_d ? flop0 : flop_d);

  if (fxt->par == NULL && mat->nrow * 2 < mat->ncol * 3) {
    /* nearly square */

    if (mat->ncol >= MINSIZE * 2 &&
	mat->nrow - mat->ncol / 2 >= MINSIZE) {
      /* try no-interpolating divide */

      flop_r = fxtld_make_div(fxt, mat, p, eps, th, flop_c);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0;

      if (fxt->ch0 != NULL)	/* success! */
	goto fin;
    }

  } else {
    /* long rectangle */

    if (mat->nrow - mat->ncol >= MINSIZE && mat->ncol >= MINSIZE) {
      /* try interpolation */

      flop_i = fxtld_make_fmm(fxt, mat, p, eps, th, flop_c);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0;

      if (fxt->fmm != NULL) {
	/* interpolation successful, choose IP-DIR or IP-DIV */

	/* get remaining matrix */
	fxt->bmat = get_scaled_bmat_ll(fxt);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0;

	/* making new set of points... */
	fxt_vecll_perm(p, fxt->idx);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0;

	/* get the new set of points */
	fxt->pp = fxt_vecll_clone(p, 0, mat->ncol - 1);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0;

	/* restore original ordering */
	fxt_vecll_iperm(p, fxt->idx);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0;

	/* check cost of IP-DIR */
	flop_id = fxt_matll_dropflop(fxt->bmat, th);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0;

	if (mat->ncol >= MINSIZE * 2) {
	  /* try IP-DIV */

	  flop_iv = fxtld_make_div(fxt, fxt->bmat, fxt->pp, eps, th,
				   MIN(flop_id, flop_c - flop_i));
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;

	  if (fxt->ch0 != NULL) {
	    /* IP + DIV was successful */

	    /* get A-scale vector */
	    fxt_vecld* ascale = fxtld_get_ascale(fxt);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return 0;

	    /* invert for rescaling */
	    fxt_vecld_einv(ascale);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return 0;

	    /* descale first child */
	    fxtld_scalerow(fxt->ch0, ascale);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return 0;

	    /* descale second child */
	    fxtld_scalerow(fxt->ch1, ascale);
	    if (fxt_error_raise() > FXT_ERROR_WARN)
	      return 0;

	    /* deallocate scaling vector */
	    fxt_vecld_del(ascale);

	    /* set resulting flop count */
	    flop_r = flop_iv + flop_i;
	    goto fin;
	  }
	}

	/* divide was not good */

	fxt_vecll_del(fxt->pp);
	fxt->pp = NULL;

	fxt_matll_del(fxt->bmat);
	fxt->bmat = NULL;

	if (flop_id + flop_i < flop_c) {
	  /* choose IP-DIR */
	  fxt_vecld *ascale;

	  /* get matrix */
	  fxt->de = get_scaled_bmat_ld(fxt);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;

	  /* drop small entries */
	  fxt_matld_drop(fxt->de, th);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;
	  
	  /* get A-scale vector */
	  ascale = fxtld_get_ascale(fxt);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;

	  /* invert for descale */
	  fxt_vecld_einv(ascale);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;

	  /* descale matrix */
	  fxt_matld_scalerow(fxt->de, ascale);
	  if (fxt_error_raise() > FXT_ERROR_WARN)
	    return 0;

	  /* deallocate scaling vector */
	  fxt_vecld_del(ascale);

	  /* set resulting flop count */
	  flop_r = flop_id + flop_i;
	  goto fin;
	}

	/* IP-DIV nor IP-DIR was chosen, unmake FMM */
	fxtld_unmake_fmm(fxt);
      }
    }
  }

  /* all trial fails, use direct evaluation */
  if (flop_d < flop0 || fxt->par == NULL) {

    /* convert matrix */
    fxt->de = fxt_matll_to_matld(mat);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    fxt_matld_drop(fxt->de, th);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    flop_r = flop_d;
  }

 fin:

#if SHOWFXTMAKE
  /* show the resulting structure */
  print_depth(fxt);
  printf(" (%ld-%ld)[%ld %ld]", fxt->mat->nrow, fxt->n_drop,
	 fxt->ui0, fxt->ui1);
  if (flop_r >= flop0 && fxt->par != NULL)
    printf(" Failed\n");
  else {
    if (fxt->fmm != NULL)
      printf(" IP");
    if (fxt->de != NULL)
      printf(" DIR");
    else
      printf(" DIV");
    printf(" %ld/%ld\n", flop_r, mat->nrow * mat->ncol);
  }
#endif

  return flop_r;
}

void fxtld_unmake(fxtld *fxt) {
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_unmake: null pointer\n");
    return;
  }

  if (fxt->de != NULL) {
    fxt_matld_del(fxt->de);
    fxt->de = NULL;
  }

  if (fxt->ch0 != NULL) {
    fxt_matll *t;

    t = fxt->ch0->mat;
    fxtld_unmake(fxt->ch0);
    fxt_matll_del(t);

    t = fxt->ch1->mat;
    fxtld_unmake(fxt->ch1);
    fxt_matll_del(t);
  }

  if (fxt->pp != NULL) {
    fxt_vecll_del(fxt->pp);
    fxt->pp = NULL;
  }

  if (fxt->bmat != NULL) {
    fxt_matll_del(fxt->bmat);
    fxt->bmat = NULL;
  }

  if (fxt->fmm != NULL)
    fxtld_unmake_fmm(fxt);

  free(fxt);
}
