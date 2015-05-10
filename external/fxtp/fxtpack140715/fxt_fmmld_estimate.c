#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  estimate values
*********************************************************************/

#define KMIN 0
#define INIT 0.2
#define KMAX 16
#define RATE 42.0

/*** compute expected error for k ***/
static double expected_error(int kest) {
  double err;  int k;

  /* error for k = KMIN */
  err = INIT;

  /* exponential decrease */
  for (k = KMIN; k < kest; k++)
    err /= RATE;

  /* the results */
  return err;
}


/*** estimate kest, kmax, and set min_np ***/
void fxt_fmmld_estimate_kmax(fxt_fmmld *fmm) {

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_estimate_kmax: null pointer\n");
    return;
  }

  /* compute the averate expansion order */
  for (fmm->kest = KMIN;  fmm->kest < KMAX &&
	 expected_error(fmm->kest) * fmm->cnorm > fmm->eps;
       fmm->kest ++);

  /* estimate maximum expansion order */
  fmm->kmax = (fmm->kest + 1) * 3;

  /* the minimum number of points for a region */
  fmm->min_np = fmm->kest * 2;
  if (fmm->min_np < 2)
    fmm->min_np = 2;
}


/********************************************************************/

static long estimate_mflop(fmmld_reg *reg,
			   fmmld_reg *mreg, fmmld_reg *lreg) {
  fmmld_cell *cell;
  int kest = reg->pfmm->kest;
  long flop = 0;

  /* M->M */
  if ((reg->type & FMM_RTMUL) && mreg != NULL)
    flop += kest * (kest + 1) / 2;

  /* L->L */
  if ((reg->type & FMM_RTLOC) && lreg != NULL)
    flop += kest * (kest + 1) / 2;

  if (reg->type & FMM_RTMUL)
    mreg = reg;

  if (reg->type & FMM_RTLOC)
    lreg = reg;

  /* P->M */
  if ((reg->type & FMM_RTEXP) && mreg != NULL) {
    int m = kest < reg->ns ? kest : reg->ns;

    flop += kest * reg->ns - m * (m - 1) / 2;

    mreg = NULL;
  }

  /* L->P */
  if ((reg->type & FMM_RTEVL) && lreg != NULL) {
    flop += kest * reg->nt;

    lreg = NULL;
  }

  /* transforms */
  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;

    if (cell->type == FMM_CTTRN)
      flop += kest * (kest + 1) / 2;

    else if (cell->type == FMM_CTDIR) {
      flop += reg->nt * creg->ns;

    } else if (cell->type == FMM_CTMUL) {
      flop += reg->nt * kest;

    } else if (cell->type == FMM_CTLOC) {
      int m = kest < creg->ns ? kest : creg->ns;

      flop += kest * creg->ns - m * (m - 1) / 2;
    }
  }

  if (reg->ch0 != NULL) {
    flop += estimate_mflop(reg->ch0, mreg, lreg);
    flop += estimate_mflop(reg->ch1, mreg, lreg);
  }

  return flop;
}

long fxt_fmmld_estimate_mflop(fxt_fmmld *fmm) {
  int i;
  long flop = 0;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_del: null pointer\n");
    return 0;
  }

  /* sum-up mflops */
  for (i=0; i< fmm->n_roots; i++)
    flop += estimate_mflop(&fmm->roots[i], NULL, NULL);

  return flop;
}

