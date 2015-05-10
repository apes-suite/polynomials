#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generate FMM for FXT
*********************************************************************/

/*** enumerate dropped rows ***/
static void make_drop(fxtld *fxt, fxt_matll *mat, double th) {
  long i, j, k0, k1;

  k0 = 0;  k1 = mat->nrow - 1;
  for (i=0; i< mat->nrow; i++) {

    /* search for large entry */
    for (j=0; j< mat->ncol && fabs(MENT(mat, i, j)) <= th; j++);

    if (j < mat->ncol)
      /* large entry found */
      fxt->idx->v[k0 ++] = i;

    else
      /* dropped row */
      fxt->idx->v[k1 --] = i;
  }

  /* set number of dropped rows */
  fxt->n_drop = mat->nrow - 1 - k1;
}


long fxtld_make_fmm(fxtld *fxt, fxt_matll *mat, fxt_vecll *x,
		    double eps, double th, long flop0) {
  fxt_matll *mat1;		/* dropped matrix */
  fxt_matld *mat2;		/* matrix in double */
  fxt_matld *bmat;		/* B-matrix */
  fxt_vecld *bvec;		/* B-vector */
  fxt_vecld *x2;		/* points in double */
  fxt_vecl *perm;		/* permutation for interpolation */
  double anorm, bnorm;		/* A-norm and B-norm */
  long flop_m, flop_s;		/* mflop for FMM */
  
  /* check null pointers */
  if (fxt == NULL || mat == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_make_fmm: null pointer\n");
    return 0;
  }

  /* initialize as no-FMM */
  fxt->n_drop = 0;
  fxt->idx = NULL;
  fxt->fmm = NULL;

  /* check over cost */
  if (flop0 <= 0)
    return 0;

  /* allocate index */
  fxt->idx = fxt_vecl_new(mat->nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* drop rows */
  make_drop(fxt, mat, th);

  /* allocate interpolation matrix */
  fxt->imat = fxt_matld_new(mat->nrow - fxt->n_drop - mat->ncol,
			    mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* number of additions */
  flop_s = mat->ncol;

  /* allocate source point vector for FMM */
  fxt->xs = fxt_vecld_new(mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* allocate target point vector for FMM */
  fxt->xt = fxt_vecld_new(mat->nrow - fxt->n_drop - mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* make dropped data */
  if (fxt->n_drop == 0) {
    mat1 = mat;
    perm = fxt->idx;
  } else {
    long i;

    /* make permutation vector */
    perm = fxt_vecl_new(mat->nrow - fxt->n_drop);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    for (i=0; i< perm->n; i++)
      perm->v[i] = i;

    /* permute matrix as dropping */
    fxt_matll_permrow(mat, fxt->idx);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* make dropped matrix */
    mat1 = fxt_matll_clone(mat, 0, mat->nrow - fxt->n_drop - 1,
			   0, mat->ncol - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* restore the original matrix */
    fxt_matll_ipermrow(mat, fxt->idx);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;
  }

  /* make interpolation matrix */
  fxt_matll_ipstab(fxt->imat, mat1, perm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* deallocate temporal matrix */
  if (mat1 != mat)
    fxt_matll_del(mat1);

  /* permute fxt->idx as interpolation */
  if (perm != fxt->idx) {
    fxt_vecl *t;

    /* get the undropped part */
    t = fxt_vecl_clone(fxt->idx, 0, fxt->idx->n - fxt->n_drop - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* permute it */
    fxt_vecl_perm(t, perm);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* store it in fxt->idx */
    fxt_vecl_prtset(fxt->idx, 0, fxt->idx->n - fxt->n_drop - 1, t);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    /* deallocate temporal */
    fxt_vecl_del(t);
    fxt_vecl_del(perm);
  }


  /* making xs and xt ... */
  x2 = fxt_vecll_to_vecld(x);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* permute x2 as interpolation */
  fxt_vecld_perm(x2, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* get source points */
  fxt_vecld_setprt(fxt->xs, x2, 0, fxt->xs->n - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* get target points */
  fxt_vecld_setprt(fxt->xt, x2, fxt->xs->n, x2->n - fxt->n_drop - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* deallocate temporal points */
  fxt_vecld_del(x2);


  /* making B-matrix ... */
  mat2 = fxt_matll_to_matld(mat);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* permute it as interpolation */
  fxt_matld_permrow(mat2, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* get B-matrix */
  bmat = fxt_matld_clone(mat2, 0, mat->ncol - 1, 0, mat->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* allocate B-vector */
  bvec = fxt_vecld_new(mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* compute B-vector */
  fxt_matld_norm2rows(bvec, bmat);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* scale interpolation matrix */
  fxt_matld_scalecol(fxt->imat, bvec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;


  /* computing B-norm ... */
  fxt_vecld_einv(bvec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* scale B-matrix */
  fxt_matld_scalerow(bmat, bvec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* compute B-norm */
  bnorm = fxt_matld_norm2(bmat, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* deallocate data */
  fxt_vecld_del(bvec);
  fxt_matld_del(bmat);
  fxt_matld_del(mat2);

  /* compute A-norm */
  anorm = fxtld_get_anorm(fxt, POWERPREC);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* make FMM ... */
  fxt->fmm = fxt_fmmld_new(fxt->imat, fxt->xs, fxt->xt,
			   eps / (anorm * bnorm));
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* check cost */
  flop_m = fxt_fmmld_evaluate_mflop(fxt->fmm);

  if (fxt_fmmld_warn(fxt->fmm) || flop_m + flop_s >= flop0 ||
      flop_m + flop_s >= fxt->imat->nrow * fxt->imat->ncol) {

    /* undo make */
    fxt_fmmld_del(fxt->fmm);
    fxt->fmm = NULL;

    fxt_vecld_del(fxt->xt);
    fxt->xt = NULL;

    fxt_vecld_del(fxt->xs);
    fxt->xs = NULL;

    fxt_matld_del(fxt->imat);
    fxt->imat = NULL;

    fxt_vecl_del(fxt->idx);
    fxt->idx = NULL;

    fxt->n_drop = 0;

    return 0;
  }

  /* rescale B-vector */
  mat2 = fxt_matll_to_matld(mat);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* permute it as interpolation */
  fxt_matld_permrow(mat2, fxt->idx);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* get B-matrix */
  bmat = fxt_matld_clone(mat2, 0, mat->ncol - 1, 0, mat->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* allocate B-vector */
  bvec = fxt_vecld_new(mat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* compute B-vector */
  fxt_matld_norm2rows(bvec, bmat);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* invert for descaling */
  fxt_vecld_einv(bvec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* scale FMM */
  fxt_fmmld_scalecol(fxt->fmm, bvec);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  /* deallocate data */
  fxt_vecld_del(bvec);
  fxt_matld_del(bmat);
  fxt_matld_del(mat2);

  /* allocate working vector */
  fxt->ws = fxt_vecld_new(fxt->imat->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  fxt->wt = fxt_vecld_new(fxt->imat->nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0;

  return flop_m + flop_s;
}

void fxtld_unmake_fmm(fxtld *fxt) {
  
  /* check null pointers */
  if (fxt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_unmake_fmm: null pointer\n");
    return;
  }

  if (fxt->fmm != NULL) {
    fxt_vecld_del(fxt->wt);
    fxt->wt = NULL;

    fxt_vecld_del(fxt->ws);
    fxt->ws = NULL;

    fxt_fmmld_del(fxt->fmm);
    fxt->fmm = NULL;

    fxt_vecld_del(fxt->xt);
    fxt->xt = NULL;

    fxt_vecld_del(fxt->xs);
    fxt->xs = NULL;

    fxt_matld_del(fxt->imat);
    fxt->imat = NULL;

    fxt_vecl_del(fxt->idx);
    fxt->idx = NULL;

    fxt->n_drop = 0;
  }
}
