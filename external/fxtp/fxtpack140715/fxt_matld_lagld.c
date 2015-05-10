#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_matld.h"
#include "fxt_lagld_loc.h"

/*********************************************************************
  matrix to linear algorithm graph
*********************************************************************/

#define MAX(x, y) ((x) > (y) ? (x) : (y))

void fxt_matld_lagld(fxt_lagld *lag, fxt_lagld_vec *y, long yini,
		     fxt_matld *mat, fxt_lagld_vec *x, long xini) {
  long i, j;

  /* check null pointers */
  if (lag == NULL || mat == NULL || y == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld: null pointer\n");
    return;
  }

  /* check index */
  if (yini < 0 || y->n <= yini || xini < 0 || x->n <= xini) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld: irregal index");
    return;
  }

  /* check size */
  if (y->n < yini + mat->nrow || x->n < xini + mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld: irregal range\n");
    return;
  }

  for (i=0; i< mat->nrow; i++)
    for (j=0; j< mat->ncol; j++)
      if (MENT(mat, i, j) != 0.0) {

	/* add an entry */
	fxt_lagld_add(lag, &y->v[i + yini],
		      MENT(mat, i, j), &x->v[j + xini]);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return;
      }
}


double fxt_matld_lagld_error(fxt_lagld *lag, fxt_matld *mat,
			     double prec, double th) {
  fxt_vecld *u, *uu, *v, *vv;
  double s, s0, ss;
  int c, cc;

  /* check null pointers */
  if (lag == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld_error: null pointer\n");
    return 0.0;
  }

  /* check size */
  if (lag->nrow != mat->nrow || lag->ncol != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld_error: size mismatch\n");
    return 0.0;
  }

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld_error: irregal precision\n");
    return 0.0;
  }

  if (lag->nrow == 0 || lag->ncol == 0)
    return 0.0;

  u  = fxt_vecld_new(lag->ncol);
  uu = fxt_vecld_new(lag->ncol);
  v  = fxt_vecld_new(lag->nrow);
  vv = fxt_vecld_new(lag->nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    do {

      fxt_vecld_rand(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

    } while (fxt_vecld_norm2(uu) == 0.0);

    s = s0 = 0.0;

    for (c=0; ; c++) {
      fxt_vecld_normalize(uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_setmatvec(v, mat, uu);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      if (th == 0.0)
	fxt_lagld_evl(vv, lag, uu);
      else
	fxt_lagld_thevl(vv, lag, uu, th);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_vecld_sub(vv, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_settrmatvec(u, mat, vv);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      if (th == 0.0)
	fxt_lagld_exp(uu, lag, vv);
      else
	fxt_lagld_thexp(uu, lag, vv, th);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_vecld_sub(uu, u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      s = fxt_vecld_norm2(uu);

      if (s == 0.0 || (c > MINPOWERITER && (s - s0) / s <= prec))
	break;

      s0 = s;
    }

    s = MAX(s, s0);

    if (s == 0.0 || fabs(s - ss) / MAX(s, ss) <= prec)
      break;

    ss = MAX(s, ss);
  }

  s = MAX(s, ss);

  fxt_vecld_del(vv);
  fxt_vecld_del(v);
  fxt_vecld_del(uu);
  fxt_vecld_del(u);

  return sqrt(s);
}


void fxt_matld_lagld_opt(fxt_lagld *lag, fxt_matld *mat, double eps) {
  double err, th;
  double err2;

  /* check null pointers */
  if (lag == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_matld_lagld_opt: null pointer\n");
    return;
  }

  /* set weight */
  fxt_lagld_setweight(lag);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* check error */
  err = fxt_matld_lagld_error(lag, mat, POWERPREC, 0.0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* for safety */
  err2 = fxt_matld_lagld_error(lag, mat, POWERPREC, 0.0);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  err = MAX(err, err2);

  /* check initial error */
  if (err >= eps * (1.0 - POWERPREC)) {
#if DISALLOWHIGHERROR
    if (err >= eps + 1e-16 * MAX(mat->nrow, mat->ncol)) {
      fxt_error_set(FXT_ERROR_USAGE,
		    "fxt_matld_lagld_opt: large initial error\n");
      return;
    } else
#endif
    {
      fxt_error_set(FXT_ERROR_WARN,
		    "fxt_matld_lagld_opt: large initial error\n");
      return;
    }
  }

  /* safety factor for the precision of the power method */
  eps *= 1.0 - POWERPREC;

  /* initial guess */
  th = eps;

  /* check error */
  err = fxt_matld_lagld_error(lag, mat, POWERPREC, th);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* for safety */
  err2 = fxt_matld_lagld_error(lag, mat, POWERPREC, th);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  err = MAX(err, err2);

  while (err > eps) {
    if (th * mat->nrow * mat->ncol < eps) {
#if DISALLOWHIGHERROR
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_matld_lagld_opt: precision unattained\n");
#else
      fxt_error_set(FXT_ERROR_WARN,
		    "fxt_matld_lagld_opt: precision unattained\n");
#endif
      break;
    }

    if (err * 0.81 > eps)
      th *= sqrt(eps / err);
    else
      th *= 0.9;

    /* check error */
    err = fxt_matld_lagld_error(lag, mat, POWERPREC, th);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    /* for safety, recheck error */
    err2 = fxt_matld_lagld_error(lag, mat, POWERPREC, th);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;

    err = MAX(err, err2);
  }

  /* zeroing small entries */
  fxt_lagld_thzero(lag, th);

  /* re-compute weight */
  fxt_lagld_setweight(lag);

  /* zeroing small entries */
  fxt_lagld_thzero(lag, 0.0);

  /* remove zeros */
  fxt_lagld_normalize(lag);
}
