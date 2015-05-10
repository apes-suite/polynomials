#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

#include "test_fxtpack.h"

#include "fxt_error.h"
#include "fxt_error.h"
#include "fxt_flptld.h"
#include "fxt_faltld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  testing FXTPACK routines
*********************************************************************/

double errmax;
static int warn_count = 0;

static void report(void) {
  if (fxt_error_level() > FXT_ERROR_WARN) {
    fprintf(stderr, "an error found: %s\n", fxt_error_message());
    exit(0);
  }

  if (fxt_error_level() == FXT_ERROR_WARN) {
    fprintf(stderr, "warning: %s\n", fxt_error_message());
    fxt_error_clear();
    warn_count ++;
  }
}

int main(int argc, char **argv) {
  long p, pp, n, m, mflop0, mflop1;
  double eps;
  fxt_vecl *mv;

  if (argc == 1) {
    printf("- test program for fxtpack\n");

    test_error();
    report();

    test_math_lsqrt();
    report();

    test_vecl();
    report();

    test_vecld();
    report();

    test_vecll();
    report();

    test_matld();
    report();

    test_matld_vecld();
    report();

    test_matll();
    report();

    test_matld_drop();
    report();

    test_matll_drop();
    report();

    test_matld_lowrank();
    report();

    test_matld_glq();
    report();

    test_matld_ipstab();
    report();

    test_matll_ipstab();
    report();

    test_math_gaussll();
    report();

    test_math_legendrell();
    report();

    test_fmmld();
    report();

    mflop0 = mflop1 = 0;

    p = 128;
    n = 85;
    eps = 1e-8;
    pp = (p + 1) / 2;
    for (m=0; m<= n; m++) {
      mflop1 += test_fxtld(p, n, m, eps);
      mflop0 += pp * (n - m + 1);
      report();
    }
    printf("- %ld %ld %e %ld (%ld: %f)\n", p, n, eps,
	   mflop1, mflop0, mflop0 / (double) mflop1);

    printf("- test file save and load...\n");

    fxt_flptld_preproc(p, p-1, eps, "test_fxtpack.flpt");
    report();

    test_flptld("test_fxtpack.flpt");
    report();

    mv = fxt_vecl_new(n+1);
    for (m=0; m<= n; m++)
      mv->v[m] = m;

    fxt_faltld_preproc(p, n, mv, eps, "test_fxtpack.falt");
    report();

    fxt_vecl_del(mv);

    test_faltld("test_fxtpack.falt");
    report();
    
  } else if (argc == 5 && '0' <= argv[4][0] && argv[4][0] <= '9') {
    errmax = 0.0;
    mflop0 = mflop1 = 0;

    p = atoi(argv[1]);
    n = atoi(argv[2]);
    m = atoi(argv[3]);
    eps = atof(argv[4]);

    pp = (p + 1) / 2;

    mflop1 += test_fxtld(p, n, m, eps);
    mflop0 += pp * (n - m + 1);
    report();

    printf("- TOTAL %ld %ld %ld %e %e %ld (%ld: %f)\n",
	   p, n, m, eps, errmax,
           mflop1, mflop0, (double) mflop0 / (double) mflop1);

  } else if (argc == 4 && '0' <= argv[3][0] && argv[3][0] <= '9') {
    errmax = 0.0;
    mflop0 = mflop1 = 0;

    p = atoi(argv[1]);
    n = atoi(argv[2]);
    eps = atof(argv[3]);

    pp = (p + 1) / 2;

    for (m=0; m <= n; m++) {
      mflop1 += test_fxtld(p, n, m, eps);
      mflop0 += pp * (n - m + 1);
      report();
    }

    printf("- TOTAL %ld %ld %e %e %ld (%ld: %f)\n", p, n, eps, errmax,
           mflop1, mflop0, (double) mflop0 / (double) mflop1);

  } else if (argc == 6) {
    long step;

    errmax = 0.0;
    mflop0 = mflop1 = 0;

    p = atoi(argv[1]);
    n = atoi(argv[2]);
    eps = atof(argv[3]);
    m = atoi(argv[4]);
    step = atoi(argv[5]);

    pp = (p + 1) / 2;

    for (; m <= n; m+=step) {
      mflop1 += test_fxtld(p, n, m, eps);
      mflop0 += pp * (n - m + 1);
      report();
    }

    printf("- TOTAL %ld %ld %e %e %ld (%ld: %f)\n", p, n, eps, errmax,
           mflop1, mflop0, (double) mflop0 / (double) mflop1);

  } else if (argc == 3) {
    errmax = 0.0;

    p = atoi(argv[1]);
    n = p - 1;
    eps = atof(argv[2]);

    pp = (p + 1) / 2;

    mflop1 = test_fxtld(p, n, 0, eps);
    mflop0 = pp * (n + 1);
    report();

    printf("- TOTAL %ld %e %e %ld (%ld: %f)\n", p, eps, errmax,
           mflop1, mflop0, (double) mflop0 / (double) mflop1);

  } else if (argc == 5) {

    p = atoi(argv[1]);
    n = atoi(argv[2]);
    eps = atof(argv[3]);

    mv = fxt_vecl_new(n+1);
    for (m=0; m <= n; m++)
      mv->v[m] = m;

    fxt_faltld_preproc(p, n, mv, eps, argv[4]);
    report();

    fxt_vecl_del(mv);

  } else if (argc == 4) {

    p = atoi(argv[1]);
    n = p - 1;
    eps = atof(argv[2]);

    fxt_flptld_preproc(p, n, eps, argv[3]);
    report();
    
  } else if (argc == 2) {
    char *s;

    for (s = argv[1]; *s != '\0' && *s != '.'; s++);

    if (*s == '\0') {
      fprintf(stderr, "filename unappropriate\n");
      report();
    }

    s ++;

    if (strcmp(s, "flpt") == 0) {

      test_flptld(argv[1]);
      report();

    } else if (strcmp(s, "falt") == 0) {

      test_faltld(argv[1]);
      report();

    } else
      fprintf(stderr, "filename unappropriate\n");
  }

  if (warn_count > 0)
    printf("- test finished, %d warning(s) found, no error found\n",
	   warn_count);
  else
    printf("- test finished, no error nor warning found\n");

  return 0;
}
