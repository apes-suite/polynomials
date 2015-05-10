#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  stable interpolation
*********************************************************************/

static void reorder(fxt_vecl *p, fxt_matld *c);

/*** get stable interpolation matrix ***/
fxt_matld* fxt_matld_ipstab(fxt_matld *a, fxt_vecl *p) {
  fxt_matld *lq;		/* LQ decomposition */
  fxt_matld *c;			/* interpolation matrix */

  /* check null pointers */
  if (a == NULL || p == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_ipstab: null pointer\n");
    return NULL;
  }

  /* check size */
  if (a->nrow <= a->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_ipstab: less rows than columns\n");
    return NULL;
  }

  /* allocate result */
  c = fxt_matld_new(a->nrow - a->ncol, a->ncol);
  if (c == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* clone the matrix */
  lq = fxt_matld_clone(a, 0, a->nrow - 1, 0, a->ncol - 1);
  if (lq == NULL) {
    fxt_error_raise();
    return NULL;
  }

  /* PLQ decomposition: P L Q = A */
  fxt_matld_glq(lq, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* permute original matrix A: A = P A */
  fxt_matld_permrow(a, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* copy matrix C = P1 A1*/
  fxt_matld_setprt(c, a, a->ncol, a->nrow - 1, 0, a->ncol - 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* C = P1 A1 Q^{-1} */
  fxt_matld_glq_mqinv(c, lq);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* C = P1 A1 Q^{-1} L0^{-1} = P1 A1 A0^{-1} P0 */
  fxt_matld_glq_mlinv(c, lq);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* original ordering */
  fxt_matld_ipermrow(a, p);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* C = P1^{-1} C P0^{-1} = A1 * A0^{-1} */
  reorder(p, c);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return NULL;

  /* deallocate LQ */
  fxt_matld_del(lq);

  return c;
}


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
