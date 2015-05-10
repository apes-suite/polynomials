#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_fxtld.h"
#include "fxt_fxtld_loc.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

   scale FXT
*********************************************************************/

/*** scale FXT ***/
void fxtld_scalerow(fxtld *fxt, fxt_vecld *s) {
  
  /* check null pointers */
  if (fxt == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxtld_scalerow: null pointer\n");
    return;
  }

  if (fxt->fmm != NULL) {
    fxt_vecld *ss;

    /* permute scaling vector as dropping/interpolation */
    fxt_vecld_perm(s, fxt->idx);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* scale for FMM */
    ss = fxt_vecld_clone(s, fxt->ws->n, fxt->ws->n + fxt->wt->n - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* scale FMM */
    fxt_fmmld_scalerow(fxt->fmm, ss);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* deallocate temporal scale vector */
    fxt_vecld_del(ss);

    /* scale for remaining computation */
    ss = fxt_vecld_clone(s, 0, fxt->ws->n - 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* restore original ordering */
    fxt_vecld_iperm(s, fxt->idx);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* let s be the remaining scaler */
    s = ss;
  }


  if (fxt->de != NULL) {

    /* scale matrix */
    fxt_matld_scalerow(fxt->de, s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  } else {

    /* scale first child */
    fxtld_scalerow(fxt->ch0, s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* scale second child */
    fxtld_scalerow(fxt->ch1, s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

  }


  if (fxt->fmm != NULL) {

    /* inverse scaling */
    fxt_vecld_einv(s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* cancell scaling of the remainings */
    fxt_fmmld_scalecol(fxt->fmm, s);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* deallocate s, originally named ss */
    fxt_vecld_del(s);
  }
}
