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

  make transform network for FMM
*********************************************************************/


/*** number of allocated cells ***/
static long nacell;


/*** number of used cells ***/
static long nucell;


/*** allocate specified number of cells ***/
static void allocate_cells(fxt_fmmld *fmm) {
  /* allocate cells */
  fmm->cells = (fmmld_cell*) malloc(sizeof(fmmld_cell) * nacell);
  if (fmm->cells == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_fmmld_make_network: allocation failed\n");
    return;
  }

  /* number of used cells */
  nucell = 0;
}


/*** deallocate cells ***/
static void deallocate_cells(fxt_fmmld *fmm) {
  /* deallocate cells */
  free(fmm->cells);
  fmm->cells = NULL;
}


/*** add a cell ***/
static void add_cell(fmmld_reg *trg, fmmld_reg *src, int type) {
  fmmld_cell *cell;

  /* check whether a cell is available */
  if (nacell <= nucell) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_fmmld_make_network: insufficient alloc\n");
    return;
  }

  /* get the next cell */
  cell = &trg->pfmm->cells[nucell ++];

  /* target region */
  cell->exp_reg = trg;

  /* source region */
  cell->evl_reg = src;

  /* connect evl_list */
  cell->evl_next = trg->evl_list;
  trg->evl_list = cell;

  /* connect exp_list */
  cell->exp_next = src->exp_list;
  src->exp_list = cell;

  /* type of the cell */
  cell->type = type;

  /* clear matrix */
  cell->m = NULL;
}


/*** expansion types ***/
#define FMM_CPMUL 1		/* partial multipole expansion */
#define FMM_CFMUL 2		/* multipole expansion */
#define FMM_CPLOC 4		/* partial local expansion */
#define FMM_CFLOC 8		/* local expansion */


/*** distance of region from a point ***/
static double DMIN(double p, double x0, double x1) {
  /* check the point is inside of the region */
  if (x0 <= p && p <= x1)
    return 0.0;

  /* which is nearer to p? */
  if (fabs(x0 - p) < fabs(x1 - p))
    return fabs(x0 - p);
  else
    return fabs(x1 - p);
}


/*** furthest distance of region from a point ***/
static double DMAX(double p, double x0, double x1) {
  /* which is farther from p? */
  if (fabs(x0 - p) > fabs(x1 - p))
    return fabs(x0 - p);
  else
    return fabs(x1 - p);
}


/*** compute possible expansion type ***/
static int expansion_type(fmmld_reg *trg, fmmld_reg *src, int mode) {
  double s0, s1, t0, t1;	/* point coordinates */
  int type;			/* expansion type */
  fxt_fmmld *fmm = trg->pfmm;	/* mother FMM structure */

  /* check existence of points */
  if (trg->nt == 0 || src->ns == 0)
    return 0;

  /* initialize type */
  type = 0;

  /* first source point */
  s0 = fmm->x->v[src->is];

  /* last source point */
  s1 = fmm->x->v[src->is + src->ns - 1];

  /* first target point */
  t0 = fmm->y->v[trg->it];

  /* last target point */
  t1 = fmm->y->v[trg->it + trg->nt - 1];

  /* check multipole expansions */
  if ((src->type & FMM_RTMUL) || (mode && (src->type & FMM_RTEXP))) {
    double size, dmin, dmax;

    /* size of the source region */
    size = DMAX(src->cs, s0, s1);

    /* distance of the target region from center of source */
    dmin = DMIN(src->cs, t0, t1);

    /* furthest distance of target from center of source */
    dmax = DMAX(src->cs, t0, t1);

    /* check for full multipole expansion */
    if (dmin >= size * 3.0)
      type |= FMM_CFMUL;

    /* check for partial multipole expansion */
    if (dmax >= size * 3.0)
      type |= FMM_CPMUL;
  }

  /* check for local expansions */
  if ((trg->type & FMM_RTLOC) || (mode && (trg->type & FMM_RTEVL))) {
    double size, dmin, dmax;

    /* size of the target region */
    size = DMAX(trg->ct, t0, t1);

    /* distance of the source region from center of target */
    dmin = DMIN(trg->ct, s0, s1);

    /* furthest distance of source from center of target */
    dmax = DMAX(trg->ct, s0, s1);

    /* check for full local expansion */
    if (dmin >= size * 3.0)
      type |= FMM_CFLOC;

    /* check for partial local expansion */
    if (dmax >= size * 3.0)
      type |= FMM_CPLOC;
  }

  /* return result */
  return type;
}


/*** whether target region should be divided before source ***/
static int divide_target(fmmld_reg *trg, fmmld_reg *src) {
  int type;			/* expansion type */
  int dc, sc;

  /* target cannot be divided if it has no child */
  if (trg->ch0 == NULL)
    return 0;

  /* source cannot be divided if it has no child */
  if (src->ch0 == NULL)
    return 1;

  /* get expansion type */
  type = expansion_type(trg, src, 0);

  /* I forgot the reasoning... */
  if (type & FMM_CFMUL)
    return 1;

  if (type & FMM_CFLOC)
    return 0;

  /* divide shallower one */
  if (trg->depth != src->depth)
    return (trg->depth < src->depth);

  /* initialize comparison via children */
  dc = sc = 0;

  /* target children */
  if (expansion_type(trg->ch0, src, 1) > type)
    dc ++;
  if (expansion_type(trg->ch1, src, 1) > type)
    dc ++;

  /* source children */
  if (expansion_type(trg, src->ch0, 1) > type)
    sc ++;
  if (expansion_type(trg, src->ch1, 1) > type)
    sc ++;

  /* compare target/source division */
  return (dc > sc);
}


/*** connect links between regions ***/
static void connect_links(fmmld_reg *trg, fmmld_reg *src) {
  int type;

  /* check existence of region */
  if (trg == NULL || src == NULL)
    return;

  /* check existence of points */
  if (trg->nt == 0 || src->ns == 0)
    return;

  /* type of possible expansion */
  type = expansion_type(trg, src, 0);

  /* multipole-local transform */
  if ((type & FMM_CFMUL) && (type & FMM_CFLOC))
    add_cell(trg, src, FMM_CTTRN);

  /* direct multipole evaluation */
  else if ((type & FMM_CFMUL) && (trg->type & FMM_RTEVL))
    add_cell(trg, src, FMM_CTMUL);

  /* direct local expansion */
  else if ((type & FMM_CFLOC) && (src->type & FMM_RTEXP))
    add_cell(trg, src, FMM_CTLOC);

  /* direct point-point evaluation */
  else if (trg->ch0 == NULL && src->ch0 == NULL)
    add_cell(trg, src, FMM_CTDIR);

  /* divide regions */
  else if (divide_target(trg, src)) {
    connect_links(trg->ch0, src);
    connect_links(trg->ch1, src);
  } else {
    connect_links(trg, src->ch0);
    connect_links(trg, src->ch1);
  }
}
    

/*** initialize region type for expansions ***/
static int init_regtype(fmmld_reg *reg) {
  int c;			/* expansion type of children */
  fxt_fmmld *fmm;		/* mother FMM structure */

  /* check void */
  if (reg == NULL)
    return 0;

  /* mother FMM structure */
  fmm = reg->pfmm;

  /* initialize children and get their type */
  c = init_regtype(reg->ch0) | init_regtype(reg->ch1);

  /* initialize region type */
  reg->type = 0;

  /* no child evaluates local exp, still I have target points */
  if ((c & FMM_RTLOC) == 0 && reg->nt > 0)
    reg->type |= FMM_RTEVL;

  /* no child expands multipole exp, still I have source points */
  if ((c & FMM_RTMUL) == 0 && reg->ns > 0)
    reg->type |= FMM_RTEXP;

  /* check for multipole expansion */
  if (reg->ns >= fmm->min_np && (fmm->approx_type & FMM_ATMUL)) {
    reg->type |= FMM_RTMUL;
    c |= FMM_RTMUL;
  }

  /* check for local expansion */
  if (reg->nt >= fmm->min_np && (fmm->approx_type & FMM_ATLOC)) {
    reg->type |= FMM_RTLOC;
    c |= FMM_RTLOC;
  }

  /* return summed expansion types */
  return c;
}


/*** remove unnecessary expansions ***/
static void modify_regtype(fmmld_reg *reg,
			   fmmld_reg *mreg, fmmld_reg *lreg) {
  /* check multipole expansion */
  if (reg->type & FMM_RTMUL) {
    fmmld_cell *cell;

    /* number of requests */
    int c = (mreg != NULL ? 1 : 0);

    /* check cells */
    for (cell = reg->exp_list; cell != NULL; cell = cell->exp_next)
      if (cell->type == FMM_CTTRN || cell->type == FMM_CTMUL)
	c ++;

    /* check number of requests */
    if (c == 0)
      reg->type -= FMM_RTMUL;
  }

  /* check local expansion */
  if (reg->type & FMM_RTLOC) {
    fmmld_cell *cell;

    /* number of requests */
    int c = (lreg != NULL ? 1 : 0);

    /* check cells */
    for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next)
      if (cell->type == FMM_CTTRN || cell->type == FMM_CTLOC)
	c ++;

    /* check number of requests */
    if (c == 0)
      reg->type -= FMM_RTLOC;
  }

  /* clear lists */
  reg->exp_list = reg->evl_list = NULL;

  /* children's mreg */
  if (reg->type & FMM_RTEXP)
    mreg = NULL;		/* multipole is expanded */
  else if (reg->type & FMM_RTMUL)
    mreg = reg;			/* I have multipole expansion */

  /* children's lreg */
  if (reg->type & FMM_RTEVL)
    lreg = NULL;		/* local expansion is evaluated */
  else if (reg->type & FMM_RTLOC)
    lreg = reg;			/* I have local expansion */

  /* recursive call */
  if (reg->ch0 != NULL) {
    modify_regtype(reg->ch0, mreg, lreg);
    modify_regtype(reg->ch1, mreg, lreg);
  }
}


/*** clear region type for expansions ***/
static void clear_regtype(fmmld_reg *reg) {
  /* check null */
  if (reg == NULL)
    return;

  /* clear lists */
  reg->evl_list = reg->exp_list = NULL;

  /* clear region types */
  reg->type = 0;

  /* recursive clear */
  clear_regtype(reg->ch0);
  clear_regtype(reg->ch1);
}


/*** make network for FMM ***/
void fxt_fmmld_make_network(fxt_fmmld *fmm) {
  int i, j;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_make_network: null pointer\n");
    return;
  }

  /* count number of possible cells */
  if ((fmm->approx_type & FMM_ATMUL) &&
      (fmm->approx_type & FMM_ATLOC)) {
    int d = fmm->max_depth;
    long n = fmm->n_roots;

    nacell = n * n + 9 * (n << d);
  } else {
    int d = fmm->max_depth;
    long n = fmm->n_roots;

    nacell = (n << d) * (n + 3 * d + 3);
  }

  /* allocate cells */
  allocate_cells(fmm);

  /* compute initial region type */
  for (i=0; i< fmm->n_roots; i++)
    init_regtype(&fmm->roots[i]);

  /* connect links */
  for (i=0; i< fmm->n_roots; i++)
    for (j=0; j< fmm->n_roots; j++)
      connect_links(&fmm->roots[i], &fmm->roots[j]);

  /* modify region type */
  for (i=0; i< fmm->n_roots; i++)
    modify_regtype(&fmm->roots[i], NULL, NULL);

  nucell = 0;

  /* reconnect links */
  for (i=0; i< fmm->n_roots; i++)
    for (j=0; j< fmm->n_roots; j++)
      connect_links(&fmm->roots[i], &fmm->roots[j]);

  /* remodify region type */
  for (i=0; i< fmm->n_roots; i++)
    modify_regtype(&fmm->roots[i], NULL, NULL);

  /* re-allocation */
  deallocate_cells(fmm);
  nacell = nucell;
  allocate_cells(fmm);

  /* reconnect links */
  for (i=0; i< fmm->n_roots; i++)
    for (j=0; j< fmm->n_roots; j++)
      connect_links(&fmm->roots[i], &fmm->roots[j]);
}


/*** undo networking ***/
void fxt_fmmld_unmake_network(fxt_fmmld *fmm) {
  int i;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_unmake_network: null pointer\n");
    return;
  }

  /* clear region types */
  for (i= fmm->n_roots - 1; i >= 0; i--)
    clear_regtype(&fmm->roots[i]);

  /* deallocate memory */
  deallocate_cells(fmm);
}
