#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

#define MIN(x, y)  ((x) < (y) ? (x) : (y))
#define MAX(x, y)  ((x) > (y) ? (x) : (y))

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  setup regions for FMM
*********************************************************************/

static void make_each_region(fxt_fmmld *fmm, fmmld_reg *reg,
			     fmmld_reg *par);


/*** setup regions ***/
void fxt_fmmld_make_regions(fxt_fmmld *fmm) {
  double xmin, xmax;
  int i;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_make_regions: null pointer\n");
    return;
  }

  /* allocate root regions */
  fmm->roots = (fmmld_reg*) malloc(sizeof(fmmld_reg) * fmm->n_roots);
  if (fmm->roots == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_fmmld_make_regions: allocation failed\n");
    return;
  }

  /* compute the extreme coordinates */
  xmin = MIN(fmm->x->v[0], fmm->y->v[0]);
  xmax = MAX(fmm->x->v[fmm->x->n - 1], fmm->y->v[fmm->y->n - 1]);

  /* xmax should be further than the last point */
  xmax += (xmax - xmin) * 1e-15;

  /* set xmin and xmax for each of root region */
  fmm->roots[0].xmin = xmin;
  fmm->roots[fmm->n_roots - 1].xmax = xmax;

  for (i= 1; i< fmm->n_roots; i++)
    fmm->roots[i-1].xmax = fmm->roots[i].xmin
      = (xmin * (fmm->n_roots - i) + xmax * i) / fmm->n_roots;

  /* recursively define regions */
  for (i=0; i< fmm->n_roots; i++)
    make_each_region(fmm, &fmm->roots[i], NULL);
}


static void unmake_each_region(fmmld_reg *reg);


/*** deallocate regions ***/
void fxt_fmmld_unmake_regions(fxt_fmmld *fmm) {
  int i;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_unmake_regions: null pointer\n");
    return;
  }

  /* deallocate recursively */
  for (i = fmm->n_roots - 1; i >= 0; i--)
    unmake_each_region(&fmm->roots[i]);

  /* deallocate root regions */
  free(fmm->roots);
  fmm->roots = NULL;
}


/*** setup each region ***/
static void setup_region(fxt_fmmld *fmm, fmmld_reg *reg,
			 fmmld_reg *par) {
  int i, j;

  /* mother structure */
  reg->pfmm = fmm;

  /* father region */
  reg->preg = par;

  /* set depth */
  if (par == NULL)
    reg->depth = 0;
  else
    reg->depth = par->depth + 1;

  /* search for the first source point */
  for (i = (par == NULL ? 0 : par->is);
       i < fmm->x->n && fmm->x->v[i] < reg->xmin; i++);
  reg->is = i;

  /* search for the last source point */
  for (j=i; j< fmm->x->n && fmm->x->v[j] < reg->xmax; j++);
  reg->ns = j - i;

  /* search for the first target point */
  for (i = (par == NULL ? 0 : par->it);
       i< fmm->y->n && fmm->y->v[i] < reg->xmin; i++);
  reg->it = i;

  /* search for the last target point */
  for (j=i; j< fmm->y->n && fmm->y->v[j] < reg->xmax; j++);
  reg->nt = j - i;

  /* center of source points */
  if (reg->ns > 0)
    reg->cs = 0.5 *
      (fmm->x->v[reg->is] + fmm->x->v[reg->is + reg->ns - 1]);
  else
    reg->cs = 0.5 * (reg->xmin + reg->xmax);

  /* center of target points */
  if (reg->nt > 0)
    reg->ct = 0.5 *
      (fmm->y->v[reg->it] + fmm->y->v[reg->it + reg->nt - 1]);
  else
    reg->ct = 0.5 * (reg->xmin + reg->xmax);

  /* initialize the other fields */
  reg->type = 0;
  reg->evl_list = reg->exp_list = NULL;

  reg->km = reg->kl = fmm->kmax;

  reg->mm = reg->ll = NULL;
  reg->pm = reg->lp = NULL;
  reg->cm = reg->cl = NULL;
}


/*** recursively create regions ***/
static void make_each_region(fxt_fmmld *fmm, fmmld_reg *reg,
			     fmmld_reg *par) {
  /* common setup */
  setup_region(fmm, reg, par);

  if (reg->depth < fmm->max_depth - 1) {
    /* make children */
    reg->ch0 = (fmmld_reg *) malloc(sizeof(fmmld_reg));
    reg->ch1 = (fmmld_reg *) malloc(sizeof(fmmld_reg));

    if (reg->ch0 == NULL || reg->ch1 == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_fmmld_make_regions: allocaion failed\n");
      return;
    }

    /* set the boundaries */
    reg->ch0->xmin = reg->xmin;
    reg->ch0->xmax =
      reg->ch1->xmin = 0.5 * (reg->xmin + reg->xmax);
    reg->ch1->xmax = reg->xmax;

    /* recursive call */
    make_each_region(fmm, reg->ch0, reg);
    make_each_region(fmm, reg->ch1, reg);
  } else
    /* no more children */
    reg->ch0 = reg->ch1 = NULL;
}


/*** deallocate regions ***/
static void unmake_each_region(fmmld_reg *reg) {
  /* recursively deallocate grand children */
  if (reg->ch1 != NULL)
    unmake_each_region(reg->ch1);
  if (reg->ch0 != NULL)
    unmake_each_region(reg->ch0);

  /* deallocate direct children */
  if (reg->ch1 != NULL) {
    free(reg->ch1);
    reg->ch1 = NULL;
  }

  if (reg->ch0 != NULL) {
    free(reg->ch0);
    reg->ch0 = NULL;
  }
}

