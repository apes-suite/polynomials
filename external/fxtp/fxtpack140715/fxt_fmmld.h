#ifndef __FXT_FMMLD_H_INCLUDED__
#define __FXT_FMMLD_H_INCLUDED__

#include "fxt_vecld.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  generalized one-dimensional fast multipole method
*********************************************************************/

/*** the data structure ***/
typedef struct FXT_FMMLD fxt_fmmld;

/*** create a new object ***/
fxt_fmmld* fxt_fmmld_new(fxt_matld *mat, fxt_vecld *xs, fxt_vecld *xt,
			 double prec);

/*** returns non-zero if the precision unattained ***/
int fxt_fmmld_warn(fxt_fmmld *fmm);

/*** deallocate memory ***/
void fxt_fmmld_del(fxt_fmmld *fmm);

/*** scale columns ***/
void fxt_fmmld_scalecol(fxt_fmmld *fmm, fxt_vecld *sc);

/*** scale rows ***/
void fxt_fmmld_scalerow(fxt_fmmld *fmm, fxt_vecld *sc);

/*** evaluation of FMM: v = mat * u ***/
void fxt_fmmld_evl(fxt_vecld *v, fxt_fmmld *fmm, fxt_vecld *u);

/*** expansion of FMM: v = mat^T * u ***/
void fxt_fmmld_exp(fxt_vecld *v, fxt_fmmld *fmm, fxt_vecld *u);

/*** get number of multiply operations ***/
long fxt_fmmld_evaluate_mflop(fxt_fmmld *fmm);

/*** get number of temporal variables ***/
long fxt_fmmld_evaluate_work(fxt_fmmld *fmm);

/*** make linear algorithm graph ***/
void fxt_fmmld_lagld(fxt_lagld *lag, fxt_lagld_vec *u,
		     fxt_fmmld *fmm, fxt_lagld_vec *v);

#endif
