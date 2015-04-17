#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"
#include "fxt_matld.h"
#include "fxt_vecld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  check preprocessing algorithm of FMM of FXTPACK
*********************************************************************/

static void test_pp(fxt_fmmld *fmm) {

  fxt_fmmld_make_regions(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  test_fmmld_regions(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_network(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  test_fmmld_network(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  test_fmmld_matrix(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_unmake_matrix(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  /* remake with true expansion orders */
  fxt_fmmld_make_matrix(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  test_fmmld_matrix(fmm);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  fxt_fmmld_unmake_matrix(fmm);

  fxt_fmmld_unmake_network(fmm);

  fxt_fmmld_unmake_regions(fmm);
}

void test_fmmld_preproc(fxt_fmmld *fmm) {
  int a, nmax, n, dmax, d;
  long np;

  fmm->eps = fmm->prec;

  fmm->kest = 5;
  fmm->kmax = 30;
  fmm->min_np = 5;

  for (a = 1; a <= 3; a ++) {
    fmm->approx_type = a;

    nmax = (fmm->ns + fmm->nt) / fmm->min_np / 2;
    if (nmax < 3)
      nmax = 3;

    for (n = 3; n <= nmax; n++) {
      fmm->n_roots = n;

      np = (fmm->ns + fmm->nt) / fmm->n_roots;
      for (dmax = 0; np >= fmm->min_np; dmax++)
	np /= 2;

      for (d = 0; d <= dmax; d ++) {
	fmm->max_depth = d;

	test_pp(fmm);
	if (fxt_error_level() > FXT_ERROR_WARN)
	  return;
      }
    }
  }
}
