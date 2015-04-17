#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test dropping for matld
*********************************************************************/

#define N 100

void test_matld_drop(void) {
  fxt_matld *a0, *a;
  long i, flop, c, prev_c;
  double eps, th;

  printf("- test_matld_drop...\n");

  a0 = fxt_matld_new(N, N);
  for (i=0; i< N*N; i++)
    a0->v[i] = exp(-60.0 * rand() / (1.0 + RAND_MAX));

  eps = 0.1;
  prev_c = 0;
  do {
    a = fxt_matld_clone(a0, 0, N-1, 0, N-1);

    th = fxt_matld_dropthreshould(a, eps);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    flop = fxt_matld_dropflop(a, th);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    fxt_matld_drop(a, th);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    for (i=c=0; i< N*N; i++) {
      if (fabs(a->v[i]) <= th && a->v[i] != 0.0) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "test_matld_drop: undropped\n");
	return;
      }

      if (fabs(a->v[i]) > th)
	c ++;
    }

    if (c != flop || c < prev_c) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_matld_drop: flop_count error\n");
      return;
    }

    fxt_matld_sub(a, a0);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    if (fxt_matld_normf(a) > eps * (1.0 + 1e-15)) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_matld_drop: norm error\n");
      return;
    }

    eps *= 0.1;
    prev_c = c;

    fxt_matld_del(a);

    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

  } while (flop < N * N);

  fxt_matld_del(a0);
}
