#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_faltld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing Fast Associated Legendre Transform
*********************************************************************/

void test_faltld_comp(fxt_faltld *falt) {
  long nm, i, j, m, wsize;
  fxt_vecld **uu, **vv, *w;
  double *u0, emax, err;
  int c;

  /* number of orders */
  nm = falt->mv->n;

  /* allocate vectors */
  uu = (fxt_vecld**) malloc(sizeof(fxt_vecld*) * (nm + 1));
  vv = (fxt_vecld**) malloc(sizeof(fxt_vecld*) * (nm + 1));
  if (uu == NULL || vv == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "test_faltld_comp: allocation failed\n");
    return;
  }

  for (i=0; i< nm; i++) {
    m = falt->mv->v[i];
    
    uu[i] = fxt_vecld_new(falt->n - m + 1);
    vv[i] = fxt_vecld_new(falt->p);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return;
  }

  /* check allocation size */
  wsize = 0;
  for (i=0; i< nm; i++) {
    /* check working array size */
    long ws = fxt_faltld_wsize(falt, falt->mv->v[i]);

    if (ws < 0) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_faltld_comp: negative work array size\n");
      return;
    }

    if (wsize < ws)
      wsize = ws;
  }

  /* allocate working array */
  w = fxt_vecld_new(wsize);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return;

  /* temporal array */
  u0 = (double *) malloc(sizeof(double) * (falt->n + 1));
  if (u0 == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "test_faltld_comp: allocation failed\n");
    return;
  }

  emax = 0;

  /* check unitarity */
  for (i=0; i< nm; i++) {
    m = falt->mv->v[i];

    for (c=0; c < 10; c++) {

      do {
	fxt_vecld_rand(uu[i]);
      } while (fxt_vecld_norm2(uu[i]) == 0.0);

      fxt_vecld_normalize(uu[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      for (j=0; j< uu[i]->n; j++)
	u0[j] = uu[i]->v[j];
      
      fxt_faltld_evl(vv[i], falt, falt->mv->v[i], uu[i], w);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      fxt_faltld_exp(uu[i], falt, falt->mv->v[i], vv[i], w);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      for (j=0; j< uu[i]->n; j++)
	uu[i]->v[j] -= u0[j];

      err = fxt_vecld_norm2(uu[i]);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return;

      if (emax < err)
	emax = err;
    }
  }

  if (emax > falt->prec) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_faltld: error %9.2e too large\n", emax);
    return;
  }

  /* deallocate temporal */
  free(u0);

  /* deallocate working array */
  fxt_vecld_del(w);

  /* deallocate vectors */
  for (i= nm - 1; i >= 0; i--) {
    fxt_vecld_del(vv[i]);
    fxt_vecld_del(uu[i]);
  }

  free(vv);
  free(uu);

  fxt_error_raise();
}
