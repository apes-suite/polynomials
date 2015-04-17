#ifndef __FXT_LA_MATLD_H_INCLUDED__
#define __FXT_LA_MATLD_H_INCLUDED__

#include "fxt_vecl.h"
#include "fxt_vecld.h"
#include "fxt_lagld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routines for matrix of double, indexed with long int
*********************************************************************/


/*** matrix data structure ***/
typedef struct {
  long nrow;			/* number of rows */
  long ncol;			/* number of columns */
  double *v;			/* entries */
} fxt_matld;


#ifndef MENT
/*** access an entry of a matrix; common the matll ***/
#define MENT(m, r, c) ((m)->v[(r) * (m)->ncol + (c)])
#endif


/*** create a new matrix ***/
fxt_matld* fxt_matld_new(long nrow, long ncol);

/*** shrink matrix ***/
void fxt_matld_shrink(fxt_matld *a, long nrow, long ncol);

/*** deallocate memory ***/
void fxt_matld_del(fxt_matld *mat);

/*** create a clone of the range ***/
fxt_matld* fxt_matld_clone(fxt_matld *mato, long rini, long rfin,
			   long cini, long cfin);

/*** fill zero to a submatrix */
void fxt_matld_zero(fxt_matld *mat, long rini, long rfin,
		    long cini, long cfin);

/*** set to identity matrix ***/
void fxt_matld_unit(fxt_matld *mat);

/*** matrix copy: A[r0:r1, c0:c1] = B ***/
void fxt_matld_prtset(fxt_matld *a,
		      long r0, long r1, long c0, long c1,
		      fxt_matld *b);

/*** matrix copy: A = B[r0:r1, c0:c1] ***/
void fxt_matld_setprt(fxt_matld *a, fxt_matld *b,
		      long r0, long r1, long c0, long c1);

/*** matrix copy: A[r0:r1, c0:c1] = B[r2:r3, c2:c3] ***/
void fxt_matld_prtsetprt(fxt_matld *a,
			 long r0, long r1, long c0, long c1,
			 fxt_matld *b,
			 long r2, long r3, long c2, long c3);

/*** matrix subtraction a - b ***/
void fxt_matld_sub(fxt_matld *a, fxt_matld *b);

/*** matrix-matrix multiplication c = a * b ***/
void fxt_matld_mul(fxt_matld *c, fxt_matld *a, fxt_matld *b);

/*** Frobenius norm ***/
double fxt_matld_normf(fxt_matld *mat);

/*** 2-norm by power method ***/
double fxt_matld_norm2(fxt_matld *mat, double prec);

/*** swap rows ***/
void fxt_matld_swaprow(fxt_matld *mat, long j, long k);

/*** swap columns ***/
void fxt_matld_swapcol(fxt_matld *mat, long j, long k);

/*** permute rows ***/
void fxt_matld_permrow(fxt_matld *mat, fxt_vecl *perm);

/*** permute columns ***/
void fxt_matld_permcol(fxt_matld *mat, fxt_vecl *perm);

/*** inverse permutation of the rows ***/
void fxt_matld_ipermrow(fxt_matld *mat, fxt_vecl *perm);

/*** get k-th row ***/
void fxt_matld_getrow(fxt_vecld *v, fxt_matld *mat, long k);

/*** get k-th column ***/
void fxt_matld_getcol(fxt_vecld *v, fxt_matld *mat, long k);

/*** scale rows (multiply diagonal matrix from right) ***/
void fxt_matld_scalerow(fxt_matld *mat, fxt_vecld *sc);

/*** scale columns (multiply diagonal matrix from left) ***/
void fxt_matld_scalecol(fxt_matld *mat, fxt_vecld *sc);

/*** 2-norm of the rows ***/
void fxt_matld_norm2rows(fxt_vecld *rnorm, fxt_matld *mat);

/*** 2-norm of the columns ***/
void fxt_matld_norm2cols(fxt_vecld *cnorm, fxt_matld *mat);


/*********************************************************************
  matrix-vector products
*********************************************************************/

/*** set matrix-vector product y = A * x ***/
void fxt_matld_setmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x);

/*** add matrix-vector product y += A * x ***/
void fxt_matld_addmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x);

/*** set matrix-partial vector product y = A * x[ini:fin] ***/
void fxt_matld_setmatprtvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x,
			    long ini, long fin);

/*** add matrix-partial vector product y += A * x[ini:fin] ***/
void fxt_matld_addmatprtvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x,
			    long ini, long fin);

/*** partial set matrix-vector product y[ini:fin] = A * x ***/
void fxt_matld_prtsetmatvec(fxt_vecld *y, long ini, long fin,
			    fxt_matld *a, fxt_vecld *x);

/*** partial add matrix-vector product y[ini:fin] += A * x ***/
void fxt_matld_prtaddmatvec(fxt_vecld *y, long ini, long fin,
			    fxt_matld *a, fxt_vecld *x);

/*** y[yini:yfin] = A * x[xini:xfin] ***/
void fxt_matld_prtsetmatprtvec(fxt_vecld *y, long yini, long yfin,
			       fxt_matld *a,
			       fxt_vecld *x, long xini, long xfin);

/*** y[yini:yfin] += A * x[xini:xfin] ***/
void fxt_matld_prtaddmatprtvec(fxt_vecld *y, long yini, long yfin,
			       fxt_matld *a,
			       fxt_vecld *x, long xini, long xfin);

/*** set transposed matrix-vector product y = A^t * x ***/
void fxt_matld_settrmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x);

/*** add transposed matrix-vector product y += A^t * x ***/
void fxt_matld_addtrmatvec(fxt_vecld *y, fxt_matld *a, fxt_vecld *x);

/*** y = A^t * x[ini:fin] ***/
void fxt_matld_settrmatprtvec(fxt_vecld *y, fxt_matld *a,
			      fxt_vecld *x, long ini, long fin);

/*** y += A^t * x[ini:fin] ***/
void fxt_matld_addtrmatprtvec(fxt_vecld *y, fxt_matld *a,
			      fxt_vecld *x, long ini, long fin);

/*** y[ini:fin] = A^t * x ***/
void fxt_matld_prtsettrmatvec(fxt_vecld *y, long ini, long fin,
			      fxt_matld *a, fxt_vecld *x);

/*** y[ini:fin] += A^t * x ***/
void fxt_matld_prtaddtrmatvec(fxt_vecld *y, long ini, long fin,
			      fxt_matld *a, fxt_vecld *x);

/*** y[yini:yfin] = A^t * x[xini:xfin] ***/
void fxt_matld_prtsettrmatprtvec(fxt_vecld *y, long yini, long yfin,
				 fxt_matld *a,
				 fxt_vecld *x, long xini, long xfin);

/*** y[yini:yfin] += A^t * x[xini:xfin] ***/
void fxt_matld_prtaddtrmatprtvec(fxt_vecld *y, long yini, long yfin,
				 fxt_matld *a,
				 fxt_vecld *x, long xini, long xfin);


/*********************************************************************
  decompositions, approximations, etc.
*********************************************************************/

/*** X * Y approximates A, returns error code ***/
int fxt_matld_lowrank(fxt_matld *a, double eps,
		      fxt_matld *x, fxt_matld *y);

/*** compute the threshold value for dropping ***/
double fxt_matld_dropthreshould(fxt_matld *mat, double eps);

/*** flop counts when dropping ***/
long fxt_matld_dropflop(fxt_matld *mat, double th);

/*** do dropping ***/
void fxt_matld_drop(fxt_matld *mat, double th);

/*** LQ decomposition by Givens rotation ***/
void fxt_matld_glq(fxt_matld *mat, fxt_vecl *piv);

/*** multiply a matrix by L^{-1} ***/
void fxt_matld_glq_mlinv(fxt_matld *mat, fxt_matld *lq);

/*** multiply a matrix by Q^{-1} ***/
void fxt_matld_glq_mqinv(fxt_matld *mat, fxt_matld *lq);

/*** get stable interpolation matrix ***/
fxt_matld* fxt_matld_ipstab(fxt_matld *a, fxt_vecl *p);

/*** make LAG (linear algorithm graph) ***/
void fxt_matld_lagld(fxt_lagld *lag, fxt_lagld_vec *y, long yini,
		     fxt_matld *mat, fxt_lagld_vec *x, long xini);

/*** compute difference of linear algrithm graph ***/
double fxt_matld_lagld_error(fxt_lagld *lag, fxt_matld *mat,
			     double prec, double th);

/*** remove small entries ***/
void fxt_matld_lagld_opt(fxt_lagld *lag, fxt_matld *mat, double prec);

#endif
