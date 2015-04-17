#ifndef __FXT_VECLD_H_INCLUDED__
#define __FXT_VECLD_H_INCLUDED__

#include <stdio.h>
#include "fxt_vecl.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  vector of double, indexed with long int
*********************************************************************/

/*** vector of double indexed with long int ***/
typedef struct {
  long n;			/* length */
  double *v;			/* vector */
} fxt_vecld;


/*** create a new vector ***/
fxt_vecld* fxt_vecld_new(long size);

/*** deallocate memory ***/
void fxt_vecld_del(fxt_vecld *x);

/*** get a clone ***/
fxt_vecld* fxt_vecld_clone(fxt_vecld *x, long ini, long fin);

/*** set value ***/
void fxt_vecld_setprt(fxt_vecld *x, fxt_vecld *y, long ini, long fin);

/*** set all zero ***/
void fxt_vecld_zero(fxt_vecld *x);

/*** set all one ***/
void fxt_vecld_one(fxt_vecld *x);

/*** set random values ***/
void fxt_vecld_rand(fxt_vecld *x);

/*** set k-th unit vector ***/
void fxt_vecld_unit(fxt_vecld *x, long k);

/*** element-wise negation ***/
void fxt_vecld_neg(fxt_vecld *x);

/*** element-wise inversion ***/
void fxt_vecld_einv(fxt_vecld *x);

/*** element-wise square ***/
void fxt_vecld_sq(fxt_vecld *x);

/*** element-wise square root ***/
void fxt_vecld_sqrt(fxt_vecld *x);

/*** scale a range ***/
void fxt_vecld_scale(fxt_vecld *x, long ini, long fin, double sc);

/*** add vector y to vector x ***/
void fxt_vecld_add(fxt_vecld *x, fxt_vecld *y);

/*** subtract vector y from vector x ***/
void fxt_vecld_sub(fxt_vecld *x, fxt_vecld *y);

/*** element-wise multiply vector y to vector x ***/
void fxt_vecld_emul(fxt_vecld *x, fxt_vecld *y);

/*** inner product of two vectors ***/
double fxt_vecld_dot(fxt_vecld *x, fxt_vecld *y);

/*** 2-norm ***/
double fxt_vecld_norm2(fxt_vecld *x);

/*** normalize in 2-norm ***/
void fxt_vecld_normalize(fxt_vecld *x);

/*** permute vector x according to perm ***/
void fxt_vecld_perm(fxt_vecld *x, fxt_vecl *perm);

/*** inverse permutation of x according to perm ***/
void fxt_vecld_iperm(fxt_vecld *x, fxt_vecl *perm);

/*** save to a file ***/
void fxt_vecld_save(FILE *f, fxt_vecld *x);

/*** restore a vector from a file ***/
fxt_vecld* fxt_vecld_restore(FILE *f);

#endif
