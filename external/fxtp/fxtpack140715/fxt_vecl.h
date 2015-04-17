#ifndef __FXT_VECL_H_INCLUDED__
#define __FXT_VECL_H_INCLUDED__

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  array of long int
*********************************************************************/

/*** array of long int ***/
typedef struct {
  long n;			/* the size of vector */
  long *v;			/* the array */
} fxt_vecl;

/*** create a new array ***/
fxt_vecl* fxt_vecl_new(long size);

/*** deallocate memory ***/
void fxt_vecl_del(fxt_vecl *x);

/*** create a clone of specified range ***/
fxt_vecl* fxt_vecl_clone(fxt_vecl *x, long ini, long fin);

/*** copy y into a range of x ***/
void fxt_vecl_prtset(fxt_vecl *x, long ini, long fin, fxt_vecl *y);

/*** swap two entries ***/
void fxt_vecl_swap(fxt_vecl *x, long i, long j);

/*** permute the vector ***/
void fxt_vecl_perm(fxt_vecl *x, fxt_vecl *perm);

#endif
