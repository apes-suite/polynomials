#ifndef __FXT_MATH_H_INCLUDED__
#define __FXT_MATH_H_INCLUDED__

#include "fxt_vecll.h"
#include "fxt_matll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  mathematical stuff
*********************************************************************/

/*** square root for long double ***/
long double fxt_lsqrt(long double);

/*********************************************************************
  special functions and related stuff
*********************************************************************/

/*** compute Gauss nodes and weights ***/
void fxt_gauss_vecll(long n, fxt_vecll **xp, fxt_vecll **wp);

/*** compute function matrix ***/
fxt_matll* fxt_legendre_matll(fxt_vecll *gv, long m, long nmax,
			      long odd);

#endif
