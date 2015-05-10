#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_fmmld.h"
#include "fxt_fmmld_loc.h"
#include "fxt_vecld.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test transform network of FMM of FXTPACK
*********************************************************************/

/*** check type of a region ***/
static void check_type(fmmld_reg *reg) {
  if ((reg->type & FMM_RTLOC) != 0 &&
      (reg->pfmm->approx_type & FMM_ATLOC) == 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_network: irregal local expansion\n");
    return;
  }

  if ((reg->type & FMM_RTLOC) != 0 && reg->nt <= 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_network: no points for local exp\n");
    return;
  }

  if ((reg->type & FMM_RTMUL) != 0 &&
      (reg->pfmm->approx_type & FMM_ATMUL) == 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_network: irregal multipole expansion\n");
    return;
  }

  if ((reg->type & FMM_RTMUL) != 0 && reg->ns <= 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_fmmld_network: no points for multipole exp\n");
    return;
  }
}

static void check_exp(fmmld_reg *reg, int type) {
  if (reg == NULL) {
    /* check whether all expansions are treated or not */

    if ((type & FMM_RTLOC) != 0 && (type & FMM_RTEVL) == 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: unevaluated local exp\n");
      return;
    }

    if ((type & FMM_RTMUL) != 0 && (type & FMM_RTEXP) == 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: unexpanded multipole exp\n");
      return;
    }

    return;
  }

  /* check expansions after evaluation or expansion */
  if (reg->type & FMM_RTEVL) {
    if ((type & FMM_RTEVL) != 0 && (reg->type & FMM_RTLOC) != 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: local exp after evl\n");
      return;
    }
  }

  if (reg->type & FMM_RTEXP) {
    if ((type & FMM_RTEXP) != 0 && (reg->type & FMM_RTMUL) != 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: multipole exp after exp\n");
      return;
    }
  }

  /* add type of this region */
  type |= reg->type;

  if (reg->ns == 0 && (type & FMM_RTMUL))
    type -= FMM_RTMUL;
  if (reg->nt == 0 && (type & FMM_RTLOC))
    type -= FMM_RTLOC;

  /* check children */
  check_exp(reg->ch0, type);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  check_exp(reg->ch1, type);
}

/* check whether all relations are treated in evaluation */
static void check_mul(fmmld_reg *reg, char *mat) {
  fmmld_cell *cell;

  if (reg == NULL)
    return;

  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *src = cell->evl_reg;
    long i, j, n;

    n = reg->pfmm->ns;

    for (i= reg->it; i < reg->it + reg->nt; i++)
      for (j= src->is; j < src->is + src->ns; j++) {

	/* check for unmark */
	if (mat[i * n + j] != '0') {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_fmmld_network: duplicated evaluation\n");
	  return;
	}

	/* mark it */
	mat[i * n + j] = '1';
      }
  }

  check_mul(reg->ch0, mat);
  check_mul(reg->ch1, mat);
}


/* check whether all relations are treated in expansion */
static void check_loc(fmmld_reg *reg, char *mat) {
  fmmld_cell *cell;

  if (reg == NULL)
    return;

  for (cell = reg->exp_list; cell != NULL; cell = cell->exp_next) {
    fmmld_reg *trg = cell->exp_reg;
    long i, j, n;

    n = reg->pfmm->ns;

    for (i= trg->it; i < trg->it + trg->nt; i++)
      for (j= reg->is; j < reg->is + reg->ns; j++) {

	/* check unmark */
	if (mat[i * n + j] != '0') {
	  fxt_error_set(FXT_ERROR_FXTBUG,
			"test_fmmld_network: duplicated expansion\n");
	  return;
	}

	/* mark it */
	mat[i * n + j] = '1';
      }
  }

  check_loc(reg->ch0, mat);
  check_loc(reg->ch1, mat);
}

void test_fmmld_network(fxt_fmmld *fmm) {
  int i;  char *mat;

  for (i=0; i< fmm->n_roots; i++) {
    check_type(&fmm->roots[i]);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  for (i=0; i< fmm->n_roots; i++) {
    check_exp(&fmm->roots[i], 0);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  mat = (char*) malloc(sizeof(char) * fmm->ns * fmm->nt);
  if (mat == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "test_fmmld_network: allocation failed\n");
    return;
  }

  for (i=0; i< fmm->ns * fmm->nt; i++)
    mat[i] = '0';

  for (i=0; i< fmm->n_roots; i++) {
    check_mul(&fmm->roots[i], mat);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  for (i=0; i< fmm->ns * fmm->nt; i++) {
    if (mat[i] != '1') {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: untreated evaluation\n");
      return;
    }

    mat[i] = '0';
  }

  for (i=0; i< fmm->n_roots; i++) {
    check_loc(&fmm->roots[i], mat);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  for (i=0; i< fmm->ns * fmm->nt; i++) {
    if (mat[i] != '1') {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_fmmld_network: untreated expansion\n");
      return;
    }
  }

  free(mat);
}
