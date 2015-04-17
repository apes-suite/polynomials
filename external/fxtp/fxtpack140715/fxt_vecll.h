#ifndef __FXT_VECLL_H_INCLUDED__
#define __FXT_VECLL_H_INCLUDED__

#include "fxt_vecl.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  vector of long double, indexed with long int
*********************************************************************/

/*** vector of long double indexed with long int ***/
typedef struct {
  long n;			/* length */
  long double *v;			/* vector */
} fxt_vecll;


/*** create a new vector ***/
fxt_vecll* fxt_vecll_new(long size);

/*** deallocate memory ***/
void fxt_vecll_del(fxt_vecll *x);

/*** get a clone ***/
fxt_vecll* fxt_vecll_clone(fxt_vecll *x, long ini, long fin);

/*** set all one ***/
void fxt_vecll_one(fxt_vecll *x);

/*** element-wise negation ***/
void fxt_vecll_neg(fxt_vecll *x);

/*** element-wise square ***/
void fxt_vecll_sq(fxt_vecll *x);

/*** element-wise square root ***/
void fxt_vecll_sqrt(fxt_vecll *x);

/*** scale a range ***/
void fxt_vecll_scale(fxt_vecll *x, long ini, long fin,
		     long double sc);

/*** element-wise multiply vector y to vector x ***/
void fxt_vecll_emul(fxt_vecll *x, fxt_vecll *y);

/*** inner product of two vectors ***/
long double fxt_vecll_dot(fxt_vecll *x, fxt_vecll *y);

/*** permute vector x according to perm ***/
void fxt_vecll_perm(fxt_vecll *x, fxt_vecl *perm);

/*** inverse permutation of x according to perm ***/
void fxt_vecll_iperm(fxt_vecll *x, fxt_vecl *perm);

/*** conversion from vecld ***/
fxt_vecll* fxt_vecld_to_vecll(fxt_vecld *x);

/*** conversion to vecld ***/
fxt_vecld* fxt_vecll_to_vecld(fxt_vecll *x);

#endif
