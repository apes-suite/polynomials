#ifndef __FXT_FLPTLD_H_INCLUDED__
#define __FXT_FLPTLD_H_INCLUDED__

#include "fxt_vecld.h"
#include "fxt_sarld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  data structure for Fast Legendre Polynomial Transform
*********************************************************************/

/*** the data structure ***/
typedef struct {
  long p;			/* number of points */
  long n;			/* degree of Legendre polynomial */
  double prec;			/* precision */

  fxt_vecld *x;			/* Gaussian nodes */
  fxt_vecld *w;			/* Gaussian weights */

  fxt_sarld *ear;		/* sarld for even half */
  fxt_sarld *oar;		/* sarld for odd  half */

  long wsize;			/* size of work array */
} fxt_flptld;

/*** create fast Legendre polynomial transform preprocessed file ***/
void fxt_flptld_preproc(long p, long n, double prec, char *fname);

/*** create fast Legendre polynomial transform data structure ***/
fxt_flptld* fxt_flptld_init(long p, long n, double prec);

/*** load fast Legendre polynomial transform ***/
fxt_flptld* fxt_flptld_load(char *fname);

/*** deallocate fast Legendre polynomial transform ***/
void fxt_flptld_del(fxt_flptld *flpt);

/*** size of working array ***/
long fxt_flptld_wsize(fxt_flptld *flpt);

/*** evaluate fast Legendre Polynomial transform ***/
void fxt_flptld_evl(fxt_vecld *v, fxt_flptld *flpt,
		    fxt_vecld *u, fxt_vecld *w);

/*** expand fast Legendre Polynomial transform ***/
void fxt_flptld_exp(fxt_vecld *u, fxt_flptld *flpt,
		    fxt_vecld *v, fxt_vecld *w);

#endif
