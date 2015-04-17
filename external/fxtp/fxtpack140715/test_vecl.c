#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "fxt_error.h"
#include "test_fxtpack.h"

#include "fxt_vecl.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test vector of long int
*********************************************************************/

static int differ(long *x, long *y, long len) {
  long i;

  for (i=0; i< len; i++)
    if (x[i] != y[i])
      return 1;

  return 0;
}

void test_vecl(void) {
  fxt_vecl *x, *y;
  long i, n;

  printf("- test_vecl...\n");

  /* create and delete */
  for (n = 0; n <= 10; n += 10) {
    x = fxt_vecl_new(n);
  
    if (x->n != n) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecl: error in new\n");
      return;
    }

    for (i=0; i< x->n; i++)
      x->v[i] = 0.0;

    fxt_vecl_del(x);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;
  }

  /* testing clone */
  n = 10;

  x = fxt_vecl_new(n);

  for (i=0; i< n; i++)
    x->v[i] = i;

  y = fxt_vecl_clone(x, 0, n-1);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, n)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecl: error in clone\n");
    return;
  }

  fxt_vecl_del(y);

  y = fxt_vecl_clone(x, 2, 7);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v + 2, y->v, 6)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecl: error in clone\n");
    return;
  }

  /* testing set */
  fxt_vecl_prtset(x, 0, 5, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 6)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecl: error in prtset\n");
    return;
  }

  fxt_vecl_prtset(x, 4, 9, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v + 4, y->v, 6)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecl: error in prtset\n");
    return;
  }

  fxt_vecl_del(y);

  /* testing swap */
  y = fxt_vecl_new(10);

  for (i=0; i< 10; i++)
    y->v[i] = i;

  for (i=0; i< 20; i++) {
    long j, k, tj, tk;

    j = (long) (10.0 * rand() / (1.0 + RAND_MAX));
    k = (long) (10.0 * rand() / (1.0 + RAND_MAX));

    if (0 > j || j >= 10 || 0 > k || k >= 10) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecl: irregal index?\n");
      return;
    }

    tj = y->v[j];
    tk = y->v[k];

    fxt_vecl_swap(y, j, k);
    if (fxt_error_level() > FXT_ERROR_WARN)
      return;

    if (tj != y->v[k] || tk != y->v[j]) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "test_vecl: error in swap\n");
      return;
    }
  }

  /* testing permute */
  for (i=0; i< 10; i++)
    x->v[i] = i;
  
  fxt_vecl_perm(x, y);
  if (fxt_error_level() > FXT_ERROR_WARN)
    return;

  if (differ(x->v, y->v, 10)) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "test_vecl: error in perm\n");
    return;
  }

  fxt_vecl_del(y);
  fxt_vecl_del(x);
}
