#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_math.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing sqrt for long double
*********************************************************************/

static double drand(void) {
  return rand() / (1.0 + RAND_MAX);
}

void test_math_lsqrt(void) {
  long double x, z, d, eps;
  int i;

  printf("- test_math_lsqrt...\n");

  /* machine epsilon for long double */
  eps = 1.0L;
  while (1.0L + eps > 1.0L)
    eps *= 0.5L;

  printf("- eps for long double = %e\n", (double) eps);

  for (i=0; i< 10; i++) {
    x = drand() + 1e-16L * drand();

    z = fxt_lsqrt(x);

    d = (z * z - x) / x;

    if (d < - 4.0L * eps || 4.0L * eps < d) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_lsqrt: error %e too large\n", (double) d);
      return;
    }
  }
}
