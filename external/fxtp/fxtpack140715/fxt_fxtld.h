#ifndef __FXT_FXTLD_H_INCLUDED__
#define __FXT_FXTLD_H_INCLUDED__

#include "fxt_matll.h"
#include "fxt_vecld.h"
#include "fxt_lagld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  preprocessing for Fast orthogonal function Transform
*********************************************************************/

/*** the data structure ***/
typedef struct FXTLD fxtld;

/*** create new data from long double matrix ***/
fxtld* fxtld_from_matll(fxt_matll *mat, fxt_vecll *p, double eps);

/*** deallocate data structure ***/
void fxtld_del(fxtld *fxt);

/*** scale rows of FXT ***/
void fxtld_scalerow(fxtld *fxt, fxt_vecld *s);

/*** evaluate FXT ***/
void fxtld_evl(fxt_vecld *v, fxtld *fxt, fxt_vecld *u);

/*** expand (transposed evaluate) FXT ***/
void fxtld_exp(fxt_vecld *v, fxtld *fxt, fxt_vecld *u);

/*** evaluate error of FXT ***/
double fxtld_error(fxtld *fxt, fxt_matld *mat, double prec);

/*** number of multiply operations ***/
long fxtld_evaluate_mflop(fxtld *fxt);

/*** number of temporal variables ***/
long fxtld_evaluate_work(fxtld *fxt);

/*** create linear algorithm graph ***/
fxt_lagld *fxtld_get_lagld(fxtld *fxt);

/*** set linear algorithm graph ***/
void fxtld_lagld(fxt_lagld *lag, fxt_lagld_vec *v,
		 fxtld *fxt, fxt_lagld_vec *u);

#endif
