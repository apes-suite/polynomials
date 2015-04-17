#ifndef __FXT_SARLD_H_INCLUDED__
#define __FXT_SARLD_H_INCLUDED__

#include <stdio.h>

#include "fxt_vecld.h"
#include "fxt_lagld.h"
#include "fxt_lagld_loc.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  Simple Array Representation for FXTPACK
*********************************************************************/

/** the data structure ***/
typedef struct {
  long ncol;			/* number of input variables */
  long nrow;			/* number of output variables */
  long ntmp;			/* number of temporal variables */

  long *p;			/* array of starting indeces for k/v */
  long *k;			/* row indeces */
  double *v;			/* matrix enries */
} fxt_sarld;

/*** allocate memory for a new data structure ***/
fxt_sarld *fxt_sarld_new(long nrow, long ncol, long ntmp, long nent);

/*** create new data stucture from LAGLD ***/
fxt_sarld *fxt_sarld_from_lagld(fxt_lagld *lag);

/*** deallocate memory ***/
void fxt_sarld_del(fxt_sarld *ar);

/*** scale by a scalar ***/
void fxt_sarld_scale(fxt_sarld *ar, double s);

/*** scale rows ***/
void fxt_sarld_scalerow(fxt_sarld *ar, fxt_vecld *sc);

/*** evaluate SARLD ***/
void fxt_sarld_evl(fxt_sarld *ar, double *w);

/*** expand (transposed evaluate) SARLD ***/
void fxt_sarld_exp(fxt_sarld *ar, double *w);

/*** compute the difference from a matrix ***/
double fxt_sarld_error(fxt_sarld *ar, fxt_matld *mat, fxt_vecld *sc,
		       double prec);

/*** save data to a file ***/
void fxt_sarld_save(FILE *f, fxt_sarld *ar);

/*** load data from a file ***/
fxt_sarld* fxt_sarld_restore(FILE *f);

#endif
