#ifndef __FXT_LAGLD_H_INCLUDED__
#define __FXT_LAGLD_H_INCLUDED__

#include "fxt_vecld.h"

/*********************************************************************
  Linear Algorithm Graph
*********************************************************************/

typedef struct FXT_LAGLD fxt_lagld;
typedef struct FXT_LAGLD_VAR fxt_lagld_var;

/*** vector of lagld_var ***/
typedef struct FXT_LAGLD_VEC {
  long n;			/* length */
  fxt_lagld_var *v;		/* array of variables */
} fxt_lagld_vec;


/*** allocate memory ***/
fxt_lagld* fxt_lagld_new(long ncol, long nrow, long nvar, long nent);

/*** deallocate memory for LAG ***/
void fxt_lagld_del(fxt_lagld *lag);

/*** get the input vector ***/
fxt_lagld_vec* fxt_lagld_ivec(fxt_lagld *lag);

/*** get the output vector ***/
fxt_lagld_vec* fxt_lagld_ovec(fxt_lagld *lag);

/*** get a new vector ***/
fxt_lagld_vec* fxt_lagld_newvec(fxt_lagld *lag, long n);

/*** add y = a * x to LAG ***/
void fxt_lagld_add(fxt_lagld *lag, fxt_lagld_var *y,
		   double a, fxt_lagld_var *x);

/*** normalize ***/
void fxt_lagld_normalize(fxt_lagld *lag);

/*** evaluate LAG ***/
void fxt_lagld_evl(fxt_vecld *v, fxt_lagld *lag, fxt_vecld *u);

/*** expand (transposed evaluation) LAG ***/
void fxt_lagld_exp(fxt_vecld *u, fxt_lagld *lag, fxt_vecld *v);

void fxt_lagld_setweight(fxt_lagld *lag);

void fxt_lagld_thzero(fxt_lagld *lag, double th);

/*** evaluate LAG with threshould ***/
void fxt_lagld_thevl(fxt_vecld *v, fxt_lagld *lag, fxt_vecld *u,
		     double th);

/*** expand (transposed evaluation) LAG with threshould ***/
void fxt_lagld_thexp(fxt_vecld *u, fxt_lagld *lag, fxt_vecld *v,
		     double th);

/*** evaluate number of multiplication ***/
long fxt_lagld_evaluate_mflop(fxt_lagld *lag);

#endif
