#ifndef __FXT_LA_MATLL_H_INCLUDED__
#define __FXT_LA_MATLL_H_INCLUDED__

#include "fxt_vecl.h"
#include "fxt_vecll.h"
#include "fxt_matld.h"


/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for matrix of long double, indexed with long int
*********************************************************************/


/*** matrix data structure ***/
typedef struct {
  long nrow;			/* number of rows */
  long ncol;			/* number of columns */
  long double *v;		/* entries */
} fxt_matll;


#ifndef MENT
/*** access an entry of a matrix; common with matld ***/
#define MENT(m, r, c) ((m)->v[(r) * (m)->ncol + (c)])
#endif

/*** create a new matrix ***/
fxt_matll* fxt_matll_new(long nrow, long ncol);

/*** deallocate memory ***/
void fxt_matll_del(fxt_matll *mat);

/*** create a clone of the range ***/
fxt_matll* fxt_matll_clone(fxt_matll *mato, long rini, long rfin,
			   long cini, long cfin);

/*** matrix copy: A[r0:r1, c0:c1] = B ***/
void fxt_matll_set(fxt_matll *a, long r0, long r1, long c0, long c1,
		   fxt_matll *b);

/*** matrix subtraction a - b ***/
void fxt_matll_sub(fxt_matll *a, fxt_matll *b);

/*** matrix-matrix multiplication c = a * b ***/
void fxt_matll_mul(fxt_matll *c, fxt_matll *a, fxt_matll *b);

/*** Frobenius norm ***/
long double fxt_matll_normf(fxt_matll *mat);

/*** permute rows ***/
void fxt_matll_permrow(fxt_matll *mat, fxt_vecl *perm);

/*** inverse permutation of the rows ***/
void fxt_matll_ipermrow(fxt_matll *mat, fxt_vecl *perm);

/*** get k-th row ***/
void fxt_matll_getrow(fxt_vecll *v, fxt_matll *mat, long k);

/*** get k-th column ***/
void fxt_matll_getcol(fxt_vecll *v, fxt_matll *mat, long k);

/*** scale rows (multiply diagonal matrix from right) ***/
void fxt_matll_scalerow(fxt_matll *mat, fxt_vecll *sc);

/*** conversion to double ***/
fxt_matld* fxt_matll_to_matld(fxt_matll *a);

/*** conversion from double ***/
fxt_matll* fxt_matld_to_matll(fxt_matld *a);


/*********************************************************************
   utilities
*********************************************************************/

/*** compute threshould for dropping ***/
double fxt_matll_dropthreshould(fxt_matll *mat, double eps);

/*** the number of non-zero entries after dropping ***/
long fxt_matll_dropflop(fxt_matll *mat, double th);

/*** stable interpolation matrix, long double only for residual ***/
void fxt_matll_ipstab(fxt_matld *c, fxt_matll *a, fxt_vecl *p);

#endif
