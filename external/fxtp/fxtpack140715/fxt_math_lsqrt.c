#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <math.h>

#include "fxt_error.h"
#include "fxt_math.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  routine for square root of a long double number
*********************************************************************/

/*** square root for long double ***/
long double fxt_lsqrt(long double x) {
  long double z;		/* 1 / sqrt(x) */

  /* check input */
  if (x < 0.0L) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_math_lsqrt: negative input\n");
    return 0.0L;
  }

  /* check 0.0 */
  if (x == 0.0L)
    return 0.0L;

  /* approximate value for sqrt(1/x) */
  z = (long double) sqrt(1.0 / (double) x);

  /* a single Newton iteration for sqrt(1/x) */
  z -= 0.5L * (x * z * z - 1.0L) * z;

  /* the second Newton iteration for sqrt(1/x) */
  z -= 0.5L * (x * z * z - 1.0L) * z;

  /* return sqrt(x) */
  return x * z;
}
