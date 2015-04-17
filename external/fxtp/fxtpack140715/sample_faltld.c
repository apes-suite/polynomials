#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_faltld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  a sample file for Fast Associated Legendre Transform
*********************************************************************/

int main(int argc, char **argv) {
  fxt_faltld *falt;
  fxt_vecld **u, **v, *w;
  long p, n, m;
  double emax;

  if (argc != 2) {
    fprintf(stderr, "usage: %s preprocessedFile\n", argv[1]);
    return 0;
  }

  /* load Fast Associated Legendre Transform from a file */
  falt = fxt_faltld_load(argv[1]);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "load from %s failed --- %s\n", argv[1],
	    fxt_error_message());
    return 0;
  }

  /* number of points */
  p = falt->p;

  /* maximum degree */
  n = falt->n;

  /* allocate vectors */
  u = (fxt_vecld**) malloc(sizeof(fxt_vecld*) * (n + 1));
  v = (fxt_vecld**) malloc(sizeof(fxt_vecld*) * (n + 1));

  if (u == NULL || v == NULL) {
    fprintf(stderr, "allocation failed\n");
    return 0;
  }

  for (m=0; m <= n; m++) {
    u[m] = fxt_vecld_new(n - m + 1);
    v[m] = fxt_vecld_new(p);
  }

  /* allocate working array */
  w = fxt_vecld_new(fxt_faltld_wsizemax(falt));

  /* set random unit vectors */
  for (m=0; m <= n; m++) {
    do {
      fxt_vecld_rand(u[m]);
    } while (fxt_vecld_norm2(u[m]) == 0.0);

    fxt_vecld_normalize(u[m]);
  }

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* wave to physical transform */
  for (m=0; m <= n; m++)
    fxt_faltld_evl(v[m], falt, m, u[m], w);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* physical to wave transform */
  for (m=0; m <= n; m++)
    fxt_faltld_exp(u[m], falt, m, v[m], w);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* check norm of the results */
  emax = 0.0;
  for (m=0; m <= n; m++) {
    double s = fxt_vecld_norm2(u[m]);

    if (emax < fabs(s - 1.0))
      emax = fabs(s - 1.0);
  }

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found: %s\n", fxt_error_message());
    return 0;
  }

  printf("maximum norm error = %e\n", emax);

  return 0;
}
