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

  test regioning of FMM for FXTPACK
*********************************************************************/

static void test_region(fmmld_reg *reg) {
  long i;

  if (reg->preg != NULL) {
    if (reg != reg->preg->ch0 && reg != reg->preg->ch1) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in parent region\n");
      return;
    }

    if (reg == reg->preg->ch0) {
      if (reg->xmin != reg->preg->xmin  ||
	  reg->xmax >= reg->preg->xmax) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_regions: error in xrange\n");
	return;
      }
    } else {
      if (reg->xmin <= reg->preg->xmin ||
	  reg->xmax != reg->preg->xmax) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_regions: error in xrange\n");
	return;
      }
    }

    if (reg->depth != reg->preg->depth + 1) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in depth\n");
      return;
    }

    if (reg->depth > reg->pfmm->max_depth) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: region too deep\n");
      return;
    }
  }

  for (i = reg->is; i < reg->is + reg->ns; i++) {
    if (reg->xmin > reg->pfmm->x->v[i] ||
	reg->pfmm->x->v[i] >= reg->xmax) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in source range\n");
      return;
    }
  }

  for (i = reg->it; i < reg->it + reg->nt; i++) {
    if (reg->xmin > reg->pfmm->y->v[i] ||
	reg->pfmm->y->v[i] >= reg->xmax) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in target range\n");
      return;
    }
  }

  if (reg->xmin > reg->cs || reg->cs > reg->xmax) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_regions: error in source center\n");
    return;
  }

  if (reg->xmin > reg->ct || reg->ct > reg->xmax) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_regions: error in target center\n");
    return;
  }

  if (reg->ns > 0) {
    if (fabs(fabs(reg->cs - reg->pfmm->x->v[reg->is]) -
	     fabs(reg->pfmm->x->v[reg->is + reg->ns -1] - reg->cs))
	> fabs(reg->xmax - reg->xmin) * 1e-14) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in source center\n");
      return;
    }
  }

  if (reg->nt > 0) {
    if (fabs(fabs(reg->ct - reg->pfmm->y->v[reg->it]) -
	     fabs(reg->ct - reg->pfmm->y->v[reg->it + reg->nt -1]))
	> fabs(reg->xmax - reg->xmin) * 1e-14) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in target center\n");
      return;
    }
  }
}

void test_fmmld_regions(fxt_fmmld *fmm) {
  int i;  long ns, nt;

  ns = nt = 0;
  for (i=0; i< fmm->n_roots; i++) {
    ns += fmm->roots[i].ns;
    nt += fmm->roots[i].nt;

    if (fmm->roots[i].depth != 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in root depth\n");
      return;
    }

    if (i > 0 && fmm->roots[i-1].xmax != fmm->roots[i].xmin) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_regions: error in root range\n");
      return;
    }

    test_region(&fmm->roots[i]);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  if (ns != fmm->ns || nt != fmm->nt) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_regions: error in root coverage\n");
    return;
  }
}
