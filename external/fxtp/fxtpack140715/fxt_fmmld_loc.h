#ifndef __FXT_FMMLD_LOC_H_INCLUDED__
#define __FXT_FMMLD_LOC_H_INCLUDED__

#include "fxt_fmmld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  FMM local header file
*********************************************************************/

/* expansions for each region */
#define FMM_RTMUL 1		/* create multipole expansion  */
#define FMM_RTLOC 2		/* create   local   expansion */
#define FMM_RTEXP 4		/* expand multipole expansion */
#define FMM_RTEVL 8		/* evaluate  local  expansion */

/* cell types */
#define FMM_CTDIR 0		/* direct evaluation */
#define FMM_CTMUL 1		/* evaluate multipole expansion */
#define FMM_CTLOC 2		/* expand local expansion */
#define FMM_CTTRN 3		/* multipole-local transform */
#define FMM_CTINV 6		/* invalid ... for check */

/* approximation types */
#define FMM_ATMUL 1		/* multipole expansion */
#define FMM_ATLOC 2		/*   local   expansion */


typedef struct FMMLD_REGION fmmld_reg;
typedef struct FMMLD_CELL   fmmld_cell;


/*** FMM data structure ***/
struct FXT_FMMLD {
  /* definitions */
  fxt_vecld *x;			/* source point coordinates */
  fxt_vecld *y;			/* target point coordinates */
  fxt_matld *mat;		/* original matrix */
  double prec;			/* requested absolute precision */

  /* derived values */
  long ns;			/* number of source points */
  long nt;			/* number of target points */
  double cnorm;			/* C-norm */

  /* the FMM data structure */
  fmmld_reg *roots;		/* the root regions */
  int warn;			/* too high precision required */

  /* algorithmic parameters */
  double eps;			/* precision parameter */
  int approx_type;		/* approximation types */
  int n_roots;			/* number of root regions */
  int max_depth;		/* height of the trees */

  /* estimated values */
  int kest;			/* estimated expansion order */
  int kmax;			/* maximum assumed expansion order */
  int min_np;			/* number of points for expansion */

  /* allocated memory */
  fmmld_cell *cells;		/* allocated cells */
};


/*** region data structure ***/
struct FMMLD_REGION {
  /* basic parameters */
  fxt_fmmld *pfmm;		/* mother structure */
  fmmld_reg *preg;		/* father region */
  int depth;			/* depth of the region */
  int type;			/* expansion type of the region */

  /* coordinates */
  double xmin, xmax;		/* range of the region */
  double cs, ct;		/* poles for source/target */
  long is, ns;			/* index/number of source points */
  long it, nt;			/* index/number of target points */
  
  /* expansions */
  int km, kl;			/* multipole/local expansion orders */
  fxt_matld *mm, *ll;		/* multipole/local transform matrix */
  fxt_matld *pm, *lp;		/* expansion/evaluation matrices */

  fmmld_cell *evl_list;		/* list of evaluation cells */
  fmmld_cell *exp_list;		/* list of expansion  cells */

  /* children */
  fmmld_reg *ch0, *ch1;		/* two children regions */

  /* work arrays */
  fxt_vecld *cm, *cl;		/* expansion coef vectors */

  /* temporal for generation */
  fxt_matld *mp;		/* multipole evaluation matrix */
  int errc;			/* lowrank error code */

  /* temporal for lagld conversion */
  fxt_lagld_vec *vm, *vl;	/* correspond to cm and cl */
};


/*** transform cell ***/
struct FMMLD_CELL {
  int type;			/* cell type */
  fxt_matld *m;			/* transform matrix */
  fmmld_reg *evl_reg;		/* region for evl_list */
  fmmld_reg *exp_reg;		/* region for exp_list */
  fmmld_cell *evl_next;		/* link for evl_list */
  fmmld_cell *exp_next;		/* link for exp_list */
};


/********************************************************************/

void fxt_fmmld_preproc(fxt_fmmld *fmm);
void fxt_fmmld_unpreproc(fxt_fmmld *fmm);

void fxt_fmmld_estimate_kmax(fxt_fmmld *fmm);
long fxt_fmmld_estimate_mflop(fxt_fmmld *fmm);

void fxt_fmmld_evaluate_kmax(fxt_fmmld *fmm);
double fxt_fmmld_error(fxt_fmmld *fmm, double prec);

void fxt_fmmld_make_regions(fxt_fmmld *fmm);
void fxt_fmmld_make_network(fxt_fmmld *fmm);
void fxt_fmmld_make_matrix(fxt_fmmld *fmm);

void fxt_fmmld_unmake_regions(fxt_fmmld *fmm);
void fxt_fmmld_unmake_network(fxt_fmmld *fmm);
void fxt_fmmld_unmake_matrix(fxt_fmmld *fmm);

void fxt_fmmld_make_mulexpmat(fmmld_reg *reg);
void fxt_fmmld_make_locexpmat(fmmld_reg *reg);

fxt_vecld* fxt_fmmld_get_scalemul(fmmld_reg *reg);
fxt_vecld* fxt_fmmld_get_scaleloc(fmmld_reg *reg);

#endif
