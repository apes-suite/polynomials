#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_sarld.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  compare matrix and simple array representation
*********************************************************************/

#define MAX(x, y) ((x) > (y) ? (x) : (y))

double fxt_sarld_error(fxt_sarld *ar, fxt_matld *mat, fxt_vecld *sc,
		       double prec) {
  fxt_vecld *u, *v;
  double s, s0, ss, *w;
  int c, cc;
  long i;

  /* check null pointers */
  if (ar == NULL || mat == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_error: null pointer\n");
    return 0.0;
  }

  /* check precision */
  if (prec <= 0.0 || 1.0 <= prec) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_error: irregal precision\n");
    return 0.0;
  }

  /* check size */
  if (ar->nrow != mat->nrow || ar->ncol != mat->ncol) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_error: size mismatch\n");
    return 0.0;
  }

  if (ar->nrow == 0 || ar->ncol == 0)
    return 0.0;

  if (sc != NULL && sc->n != mat->nrow) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_sarld_error: size mismatch\n");
    return 0.0;
  }

  u  = fxt_vecld_new(ar->ncol);
  v  = fxt_vecld_new(ar->nrow);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return 0.0;

  w = (double *) malloc(sizeof(double) * 
			(ar->ncol + ar->nrow + ar->ntmp));
  if (w == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_sarld_error: allocation failed\n");
    return 0.0;
  }

  ss = 0.0;

  for (cc=0; cc< MAXPOWERTRY; cc++) {

    do {

      fxt_vecld_rand(u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

    } while (fxt_vecld_norm2(u) == 0.0);

    s = s0 = 0.0;

    for (c=0; ; c++) {

      fxt_vecld_normalize(u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      fxt_matld_setmatvec(v, mat, u);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      for (i=0; i< ar->ncol; i++)
	w[i] = u->v[i];

      fxt_sarld_evl(ar, w);

      for (i=0; i< ar->nrow; i++)
	v->v[i] -= w[i + ar->ncol];

      if (sc != NULL) {
	fxt_vecld_emul(v, sc);
	if (fxt_error_raise() > FXT_ERROR_WARN)
	  return 0.0;
      }

      fxt_matld_settrmatvec(u, mat, v);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return 0.0;

      for (i=0; i< ar->nrow; i++)
	w[i + ar->ncol] = v->v[i];

      fxt_sarld_exp(ar, w);

      for (i=0; i< ar->ncol; i++)
	u->v[i] -= w[i];

      s = fxt_vecld_norm2(u);

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

  free(w);
  fxt_vecld_del(v);
  fxt_vecld_del(u);

  return sqrt(s);
}
