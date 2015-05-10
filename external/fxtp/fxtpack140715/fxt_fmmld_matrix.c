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

  make matrix for FMM
*********************************************************************/

/*** allocate matrices ***/
static void alloc_mat(fmmld_reg *reg,
		      fmmld_reg *mreg, fmmld_reg *lreg) {
  fmmld_cell *cell;

  /* clear error code */
  reg->errc = 0;

  /* for multipole expansion */
  if (reg->type & FMM_RTMUL) {
    if (mreg != NULL) {
      /* multipole-multipole transform */
      reg->mm = fxt_matld_new(mreg->km, reg->km);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    } else
      reg->mm = NULL;

    /* I have multipole expansion */
    mreg = reg;
  }

  /* point-to-multipole expansion done */
  if (mreg != NULL && (reg->type & FMM_RTEXP)) {
    /* from point to multipole expansion */
    reg->pm = fxt_matld_new(reg->km, reg->ns);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* no more expansion */
    mreg = NULL;
  } else
    reg->pm = NULL;

  /* for local expansion */
  if (reg->type & FMM_RTLOC) {
    if (lreg != NULL) {
      /* local-local transform */
      reg->ll = fxt_matld_new(reg->kl, lreg->kl);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    } else
      reg->ll = NULL;

    /* I have local expansion */
    lreg = reg;
  }

  /* local-to-point expansion done */
  if (lreg != NULL && (reg->type & FMM_RTEVL)) {
    /* local to point evaluation */
    reg->lp = fxt_matld_new(reg->nt, lreg->kl);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* no more expansion */
    lreg = NULL;
  } else
    reg->lp = NULL;

  /* for transforms */
  for (cell = reg->evl_list; cell != NULL; cell = cell->evl_next) {
    fmmld_reg *creg = cell->evl_reg;
    
    if (cell->type == FMM_CTTRN) /* multipole to local transform */
      cell->m = fxt_matld_new(reg->kl, creg->km);
    else if (cell->type == FMM_CTLOC) /* direct local expansion */
      cell->m = fxt_matld_new(reg->kl, creg->ns);
    else if (cell->type == FMM_CTDIR) /* point-to-point evaluation */
      cell->m = fxt_matld_new(reg->nt, creg->ns);
    else if (cell->type == FMM_CTMUL) /* direct multipole evaluate */
      cell->m = fxt_matld_new(reg->nt, creg->km);
    else
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_fmmld_make_matrix: network error\n");

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* recursive call to children */
  if (reg->ch0 != NULL) {
    alloc_mat(reg->ch0, mreg, lreg);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
    alloc_mat(reg->ch1, mreg, lreg);
  }
}


/*** deallocate cell's matrices ***/
static void dealloc_cellmat(fmmld_cell *cell) {
  if (cell == NULL)
    return;

  /* last allocated, first deallocate */
  dealloc_cellmat(cell->evl_next);

  /* deallocate the matrix */
  fxt_matld_del(cell->m);
  cell->m = NULL;
}


/*** deallocate matrices ***/
static void dealloc_mat(fmmld_reg *reg) {
  /* recursion */
  if (reg->ch0 != NULL) {
    dealloc_mat(reg->ch1);
    dealloc_mat(reg->ch0);
  }

  /* deallocate matrices in the evl_list */
  dealloc_cellmat(reg->evl_list);

  /* simply deallocate if exists */
  if (reg->lp != NULL) {
    fxt_matld_del(reg->lp);
    reg->lp = NULL;
  }

  if (reg->ll != NULL) {
    fxt_matld_del(reg->ll);
    reg->ll = NULL;
  }

  if (reg->pm != NULL) {
    fxt_matld_del(reg->pm);
    reg->pm = NULL;
  }

  if (reg->mm != NULL) {
    fxt_matld_del(reg->mm);
    reg->mm = NULL;
  }
}


/*** allocate temporal matrices ***/
static void alloc_temp(fmmld_reg *reg) {
  if (reg->type & FMM_RTMUL) {
    /* multipole 'evaluation' matrix */
    reg->mp = fxt_matld_new(reg->pfmm->nt, reg->km);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  } else
    reg->mp = NULL;

  if ((reg->type & FMM_RTLOC) && (reg->type & FMM_RTEVL) == 0) {
    /* local evaluation matrix for non-evaluating region */
    reg->lp = fxt_matld_new(reg->nt, reg->kl);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* recursion to children */
  if (reg->ch0 != NULL) {
    alloc_temp(reg->ch0);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
    alloc_temp(reg->ch1);
  }
}


/*** deallocate temporal matrices ***/
static void dealloc_temp(fmmld_reg *reg) {
  /* recursion to children */
  if (reg->ch0 != NULL) {
    dealloc_temp(reg->ch1);
    dealloc_temp(reg->ch0);
  }

  /* remove local evaluation matrix */
  if ((reg->type & FMM_RTLOC) && (reg->type & FMM_RTEVL) == 0) {
    fxt_matld_del(reg->lp);
    reg->lp = NULL;
  }

  /* remove multipole evaluation matrix */
  if (reg->type & FMM_RTMUL) {
    fxt_matld_del(reg->mp);
    reg->mp = NULL;
  }
}


/*** allocate expansion vectors ***/
static void alloc_vec(fmmld_reg *reg) {
  if (reg->type & FMM_RTMUL) {
    /* allocate multipole expansion vector */
    reg->cm = fxt_vecld_new(reg->km);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  } else
    reg->cm = NULL;

  if (reg->type & FMM_RTLOC) {
    /* allocate local expansion vector */
    reg->cl = fxt_vecld_new(reg->kl);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  } else
    reg->cl = NULL;

  /* recursion to children */
  if (reg->ch0 != NULL) {
    alloc_vec(reg->ch0);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
    alloc_vec(reg->ch1);
  }
}


/*** deallocate expansion vectors ***/
static void dealloc_vec(fmmld_reg *reg) {
  /* recursion to children */
  if (reg->ch0 != NULL) {
    dealloc_vec(reg->ch1);
    dealloc_vec(reg->ch0);
  }

  /* deallocate local expansion vector */
  if (reg->cl != NULL) {
    fxt_vecld_del(reg->cl);
    reg->cl = NULL;
  }

  /* deallocate multipole expansion vector */
  if (reg->cm != NULL) {
    fxt_vecld_del(reg->cm);
    reg->cm = NULL;
  }
}


/*** check error codes ***/
static int check_errc(fmmld_reg *reg) {
  int errc = reg->errc;

  reg->errc = 0;

  if (errc & 1)
    reg->km ++;
  if (errc & 2)
    reg->kl ++;

  if (reg->ch0 != NULL) {
    errc |= check_errc(reg->ch0);
    errc |= check_errc(reg->ch1);
  }

  return errc;
}


/*** generate matrices ***/
void fxt_fmmld_make_matrix(fxt_fmmld *fmm) {

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_make_matrix: null pointer\n");
    return;
  }

  while (1) {
    int i, errc;

    /* allocate expansion/transform/evaluation matrices */
    for (i=0; i< fmm->n_roots; i++) {
      alloc_mat(&fmm->roots[i], NULL, NULL);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }

    /* allocate temporal matrices */
    for (i=0; i< fmm->n_roots; i++) {
      alloc_temp(&fmm->roots[i]);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }

    /* make multipole expansions */
    for (i=0; i< fmm->n_roots; i++) {
      fxt_fmmld_make_mulexpmat(&fmm->roots[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    /* make local expansions */
    for (i=0; i< fmm->n_roots; i++) {
      fxt_fmmld_make_locexpmat(&fmm->roots[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;
    }

    /* deallocate temporal matrices */
    for (i = fmm->n_roots - 1; i >= 0; i--)
      dealloc_temp(&fmm->roots[i]);

    /* allocate expansion vectors */
    for (i=0; i< fmm->n_roots; i++) {
      alloc_vec(&fmm->roots[i]);
      if (fxt_error_level() > FXT_ERROR_WARN)
	return;
    }

    /* check error code */
    errc = 0;
    for (i=0; i< fmm->n_roots; i++)
      errc |= check_errc(&fmm->roots[i]);

    if (errc == 0)
      break;

    /* retry */
    fxt_fmmld_unmake_matrix(fmm);
  }
}


/*** deallocate matrices ***/
void fxt_fmmld_unmake_matrix(fxt_fmmld *fmm) {
  int i;

  /* check null pointers */
  if (fmm == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_fmmld_unmake_matrix: null pointer\n");
    return;
  }

  /* deallocate expansion vectors */
  for (i = fmm->n_roots - 1; i >= 0; i--)
    dealloc_vec(&fmm->roots[i]);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  /* deallocate matrices */
  for (i = fmm->n_roots - 1; i >= 0; i--)
    dealloc_mat(&fmm->roots[i]);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;
}
