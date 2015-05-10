#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_flptld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  a sample file for Fast Legendre Polynomial Transform
*********************************************************************/

int main(int argc, char **argv) {
  fxt_flptld *flpt;
  fxt_vecld *u, *v, *w;
  long p, n;
  double err;

  if (argc != 2) {
    fprintf(stderr, "usage: %s preprocessedFile\n", argv[0]);
    return 0;
  }

  /* load Fast Legendre Polynomial Transform from a file */
  flpt = fxt_flptld_load(argv[1]);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "load from %s failed --- %s\n", argv[1],
	    fxt_error_message());
    return 0;
  }

  /* number of points */
  p = flpt->p;

  /* maximum degree */
  n = flpt->n;

  /* allocate vectors */
  u = fxt_vecld_new(n + 1);
  v = fxt_vecld_new(p);

  /* allocate working array */
  w = fxt_vecld_new(fxt_flptld_wsize(flpt));

  /* set random unit vector */
  do {
    fxt_vecld_rand(u);
  } while (fxt_vecld_norm2(u) == 0.0);
  
  fxt_vecld_normalize(u);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* wave to physical transform */
  fxt_flptld_evl(v, flpt, u, w);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* physical to wave transform */
  fxt_flptld_exp(u, flpt, v, w);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found --- %s\n", fxt_error_message());
    return 0;
  }

  /* check norm of the results */
  err = fabs(fxt_vecld_norm2(u) - 1.0);

  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "error found: %s\n", fxt_error_message());
    return 0;
  }

  printf("norm error = %e\n", err);

  return 0;
}
