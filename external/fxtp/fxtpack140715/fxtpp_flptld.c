#include <stdio.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "fxt_flptld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  create preprocessed file for Fast Legendre Polynomial Transform
*********************************************************************/

int main(int argc, char **argv) {
  long p, n;  double prec;

  if (argc != 5) {
    fprintf(stderr, "usage: %s #points maxdeg prec file\n",
	    argv[0]);
    exit(1);
  }

  p = atoi(argv[1]);
  n = atoi(argv[2]);
  prec = atof(argv[3]);

  if (p < 0) {
    fprintf(stderr, "negative number of points %ld\n", p);
    exit(1);
  }

  if (n < 0) {
    fprintf(stderr, "negative maximum degree %ld\n", n);
    exit(1);
  }

  if (p <= n) {
    fprintf(stderr, "more points than maximum degree are needed\n");
    exit(1);
  }

  if (prec >= 1.0) {
    fprintf(stderr, "precision %e too low\n", prec);
    exit(1);
  }

  fxt_flptld_preproc(p, n, prec, argv[4]);
  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "preprocessing failed with error: %s\n",
	    fxt_error_message());
    exit(1);
  }

  return 0;
}
