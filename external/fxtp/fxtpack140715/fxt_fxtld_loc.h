#ifndef __FXT_FXTLD_LOC_H_INCLUDED__
#define __FXT_FXTLD_LOC_H_INCLUDED__

#include "fxt_fxtld.h"
#include "fxt_vecl.h"
#include "fxt_vecld.h"
#include "fxt_matld.h"
#include "fxt_matll.h"
#include "fxt_fmmld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  internal header file for Fast orthogonal function Transform
*********************************************************************/

/*** the recursive data structure ***/
struct FXTLD {
  fxtld *par;			/* parent fxt */
  long ui0, ui1;		/* range of u-vector */
  fxt_matll *mat;		/* original matrix */
  fxt_vecll *p;			/* original point vector */

  long n_drop;			/* number of dropped rows */
  fxt_vecl *idx;		/* reordering index */

  fxt_fmmld *fmm;		/* interpolation matrix */
  fxt_matld *imat;		/* save matrix for FMM */
  fxt_vecld *xs, *xt;		/* saved points for FMM */

  fxt_matld *de;		/* direct evaluation */

  fxtld *ch0, *ch1;		/* children */
  fxt_matll *bmat;		/* saved matrix for DIV */
  fxt_vecll *pp;		/* saved point for DIV */

  fxt_vecld *ws, *wt;		/* working vector */
};

long fxtld_make(fxtld *fxt, fxt_matll *mat, fxt_vecll *p,
		double eps, double th, long flop0, long ui0);

void fxtld_unmake(fxtld *fxt);

long fxtld_make_fmm(fxtld *fxt, fxt_matll *mat, fxt_vecll *p,
		    double eps, double th, long flop0);

void fxtld_unmake_fmm(fxtld *fxt);

long fxtld_make_div(fxtld *fxt, fxt_matll *mat, fxt_vecll *p,
		    double eps, double th, long flop0);

/*** get A-norm ***/
double fxtld_get_anorm(fxtld *fxt, double prec);

/*** get A-scale ***/
fxt_vecld *fxtld_get_ascale(fxtld *fxt);

#endif
