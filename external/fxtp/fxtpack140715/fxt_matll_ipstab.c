#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matll.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  stable interpolation for long double matrix
*********************************************************************/

static void reorder(fxt_vecl *p, fxt_matld *c);

/*** stable interpolation matrix, long double only for residual ***/
void fxt_matll_ipstab(fxt_matld *c, fxt_matll *a, fxt_vecl *p) {
  fxt_matld *lq;		/* LQ decomposition in double */
  fxt_matld *r;			/* residual/update in double */
  fxt_matll *a0, *a1;		/* upper/lower parts of A */
  fxt_matll *cc, *rr;		/* long double versions of c and r */

  /* check null pointers */
  if (c == NULL || a == NULL || p == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_ipstab: null pointer\n");
    return;
  }

  /* check size */
  if (a->nrow <= a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_ipstab: rows fewer than columns\n");
    return;
  }

  /* check size of c */
  if (c->nrow != a->nrow - a->ncol || c->ncol != a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matll_ipstab: size mismatch\n");
    return;
  }

  /* make clone in double */
  lq = fxt_matll_to_matld(a);
  if (lq == NULL) {
    fxt_error_raise();
    return;
  }

  /* PLQ decomposition */
  fxt_matld_glq(lq, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* C = L1 (lower part of L) */
  fxt_matld_setprt(c, lq, a->ncol, a->nrow - 1, 0, a->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* C = L1 L0^{-1} */
  fxt_matld_glq_mlinv(c, lq);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* permute original matrix as P */
  fxt_matll_permrow(a, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* upper part of the original matrix */
  a0 = fxt_matll_clone(a, 0, a->ncol - 1, 0, a->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN) {
    return;
  }

  /* lower part of the original matrix */
  a1 = fxt_matll_clone(a, a->ncol, a->nrow - 1, 0, a->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* restore original ordering */
  fxt_matll_ipermrow(a, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* raise interpolation matrix to long double */
  cc = fxt_matld_to_matll(c);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* allocate residual R = C A0 - A1 */
  rr = fxt_matll_new(a->nrow - a->ncol, a->ncol);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* R = C A0 */
  fxt_matll_mul(rr, cc, a0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* R = C A0 - A1 */
  fxt_matll_sub(rr, a1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* lower the residual to double */
  r = fxt_matll_to_matld(rr);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* R = (C A0 - A1) Q^{-1} */
  fxt_matld_glq_mqinv(r, lq);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* R = (C A0 - A1) A0^{-1} */
  fxt_matld_glq_mlinv(r, lq);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* C = C - (C A0 - A1) A0^{-1} */
  fxt_matld_sub(c, r);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* deallocate memory */
  fxt_matld_del(r);
  fxt_matll_del(rr);
  fxt_matll_del(cc);
  fxt_matll_del(a1);
  fxt_matll_del(a0);
  fxt_matld_del(lq);

  /* reorder */
  reorder(p, c);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;
}


/*** sort each part of p, reorder c accordingly ***/
static void reorder(fxt_vecl *p, fxt_matld *c) {
  long *pos;			/* inverse of p */
  fxt_vecl *split;		/* split reordering vector */
  fxt_vecl *colp, *rowp;	/* col/row reordering vectors */
  long i, j, k;

  /* check size */
  if (p->n != c->nrow + c->ncol) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_ipstab: size mismatch in reordering\n");
    return;
  }

  /* allocate inverse index */
  pos = (long *) malloc(sizeof(long) * p->n);
  if (pos == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_matld_ipstab: allocation failed\n");
    return;
  }

  /* make inverse index */
  for (i=0; i< p->n; i++)
    pos[p->v[i]] = i;

  /* allocate split ordering vector */
  split = fxt_vecl_new(p->n);
  if (split == NULL) {
    fxt_error_raise();
    return;
  }

  /* split orders */
  j = 0;  k = c->ncol;
  for (i=0; i< p->n; i++)
    if (pos[i] < c->ncol)
      split->v[j++] = pos[i];
    else
      split->v[k++] = pos[i];

  /* check result */
  if (j != c->ncol || k != p->n) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_matld_ipstab: error in split\n");
    return;
  }

  /* reorder p */
  fxt_vecl_perm(p, split);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* check reordering */
  for (i=1; i< p->n; i++)
    if (i != c->ncol && p->v[i-1] >= p->v[i]) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_matld_ipstab: error in reordering\n");
      return;
    }

  /* column reordering vector */
  colp = fxt_vecl_clone(split, 0, c->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* reorder columns */
  fxt_matld_permcol(c, colp);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* row reordering vector */
  rowp = fxt_vecl_clone(split, c->ncol, c->ncol + c->nrow - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  for (i=0; i< c->nrow; i++)
    rowp->v[i] -= c->ncol;

  /* reorder rows */
  fxt_matld_permrow(c, rowp);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* deallocate temporals */
  fxt_vecl_del(rowp);
  fxt_vecl_del(colp);
  fxt_vecl_del(split);
  free(pos);
}
