#ifndef __FXT_FALTLD_H_INCLUDED__
#define __FXT_FALTLD_H_INCLUDED__

#include "fxt_vecl.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  Fast Associated Legendre Transform for spherical harmonics
*********************************************************************/

typedef struct FALTLD_ALT faltld_alt;

/*** the data structure ***/
typedef struct {
  long p;			/* number of points */
  long n;			/* maximum degree */
  fxt_vecl *mv;			/* array of order m */
  double prec;			/* precision */

  fxt_vecld *x;			/* Gaussian nodes */
  fxt_vecld *w;			/* Gaussian weigths */

  faltld_alt **alt;		/* array of Legendre transforms */
} fxt_faltld;

/*** create fast spherical harmonic transform preprocessed file ***/
void fxt_faltld_preproc(long p, long n, fxt_vecl *mv, double prec, 
			char *fname);

/*** load fast spherical harmonic transform ***/
fxt_faltld* fxt_faltld_load(char *fname);

/*** deallocate fast spherical harmonic transform ***/
void fxt_faltld_del(fxt_faltld *falt);

/*** size of working array ***/
long fxt_faltld_wsize(fxt_faltld *falt, long m);

/*** maximum size of working array ***/
long fxt_faltld_wsizemax(fxt_faltld *falt);

/*** evaluate fast spherical harmonic transform ***/
void fxt_faltld_evl(fxt_vecld *v, fxt_faltld *falt, long m,
		    fxt_vecld *u, fxt_vecld *w);

/*** expand fast spherical harmonic transform ***/
void fxt_faltld_exp(fxt_vecld *u, fxt_faltld *falt, long m,
		    fxt_vecld *v, fxt_vecld *w);

#endif
