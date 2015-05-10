#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_matll.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test dropping for matll
*********************************************************************/

#define N 100

void test_matll_drop(void) {
  fxt_matll *a0, *a;
  long i, flop, prev_flop;
  double eps, th;

  printf("- test_matll_drop...\n");

  a0 = fxt_matll_new(N, N);
  for (i=0; i< N * N; i++)
    a0->v[i] = (long double) exp(-60.0 * rand() / (1.0 + RAND_MAX));

  eps = 0.1;
  prev_flop = 0;
  do {
    a = fxt_matll_clone(a0, 0, N-1, 0, N-1);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    th = fxt_matll_dropthreshould(a, eps);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    flop = fxt_matll_dropflop(a, th);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    if (flop < prev_flop) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_matll_drop: abnormal number of nonzero\n");
      return;
    }

    fxt_matll_del(a);

    eps *= 0.1;
    prev_flop = flop;
  } while (flop < N * N);

  fxt_matll_del(a0);
}
