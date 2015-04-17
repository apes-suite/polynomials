#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing estimation routines for FMM of FXTPACK
*********************************************************************/

void test_fmmld_estimate(fxt_fmmld *fmm) {
  int d;  long n;

  fmm->eps = fmm->prec;

  /* do estimation */
  fxt_fmmld_estimate_kmax(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fmm->approx_type = 3;
  fmm->n_roots = 3;

  n = ((fmm->ns > fmm->nt) ? fmm->ns : fmm->nt) / 3;
  for (d = 0; n >= fmm->min_np; d ++, n /= 2);

  fmm->max_depth = d;

  /* just check to run */

  fxt_fmmld_make_regions(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_network(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_evaluate_kmax(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_unmake_matrix(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_unmake_network(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_unmake_regions(fmm);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}
