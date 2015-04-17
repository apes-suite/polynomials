#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "test_fxtpack.h"

#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test matrix generation for FMM of FXTPACK
*********************************************************************/

static void check_matrix_size(fmmld_reg *reg,
			      fmmld_reg *mreg, fmmld_reg *lreg) {
  fmmld_cell *cell;

  if ((reg->type & FMM_RTMUL) && mreg != NULL) {
    if (reg->mm == NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: MM non-existent\n");
      return;
    }

    if (reg->mm->nrow != mreg->km) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: MM row size mismatch\n");
      return;
    }

    if (reg->mm->ncol != reg->km) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: MM column size mismatch\n");
      return;
    }
  }

  if (reg->type & FMM_RTMUL)
    mreg = reg;

  if ((reg->type & FMM_RTEXP) && mreg != NULL) {
    if (reg->pm == NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: PM non-existent\n");
      return;
    }

    if (reg->pm->nrow != mreg->km) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: PM row size mismatch\n");
      return;
    }

    if (reg->pm->ncol != reg->ns) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: PM column size mismatch\n");
      return;
    }

    mreg = NULL;
  }


  if ((reg->type & FMM_RTLOC) && lreg != NULL) {
    if (reg->ll == NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LL non-existent\n");
      return;
    }

    if (reg->ll->nrow != reg->kl) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LL row size mismatch\n");
      return;
    }

    if (reg->ll->ncol != lreg->kl) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LL column size mismatch\n");
      return;
    }
  }

  if (reg->type & FMM_RTLOC)
    lreg = reg;

  if ((reg->type & FMM_RTEVL) && lreg != NULL) {
    if (reg->lp == NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LP non-existent\n");
      return;
    }

    if (reg->lp->nrow != reg->nt) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LP row size mismatch\n");
      return;
    }

    if (reg->lp->ncol != lreg->kl) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: LP column size mismatch\n");
      return;
    }

    lreg = NULL;
  }


  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;

    if (cell->m == NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_matrix: cell matrix non-existent\n");
      return;
    }

    if (cell->type == FMM_CTTRN) {
      if (cell->m->nrow != reg->kl) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: ML row size mismatch\n");
	return;
      }

      if (cell->m->ncol != creg->km) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: ML column size mismatch\n");
	return;
      }

    } else if (cell->type == FMM_CTLOC) {
      if (cell->m->nrow != reg->kl) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: PL row size mismatch\n");
	return;
      }

      if (cell->m->ncol != creg->ns) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: PL column size mismatch\n");
	return;
      }
   
    } else if (cell->type == FMM_CTDIR) {
      if (cell->m->nrow != reg->nt) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: PP row size mismatch\n");
	return;
      }

      if (cell->m->ncol != creg->ns) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: PP column size mismatch\n");
	return;
      }

    } else if (cell->type == FMM_CTMUL) {
      if (cell->m->nrow != reg->nt) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: MP row size mismatch\n");
	return;
      }

      if (cell->m->ncol != creg->km) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_fmmld_matrix: MP column size mismatch\n");
	return;
      }
    }
  }

  if (reg->ch0 != NULL) {
    check_matrix_size(reg->ch0, mreg, lreg);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    check_matrix_size(reg->ch1, mreg, lreg);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }
}

void test_fmmld_matrix(fxt_fmmld *fmm) {
  int i;
  double s;

  for (i=0; i< fmm->n_roots; i++) {
    check_matrix_size(&fmm->roots[i], NULL, NULL);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  s = fxt_fmmld_error(fmm, POWERPREC);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (s > fmm->prec * 20.0 &&
      s > 1e-16 * fmm->cnorm * (fmm->ns + fmm->nt)) {
    fxt_error_set(FXT_ERROR_WARN,
		  "test_fmmld_matrix: error %e seems large\n", s);
  }
}
