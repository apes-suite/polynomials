#ifndef __TEST_FXTPACK_H_INCLUDED__
#define __TEST_FXTPACK_H_INCLUDED__

#include "fxt_fmmld.h"
#include "fxt_faltld.h"
#include "fxt_sarld.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  test routine of fxtpack
*********************************************************************/

extern double errmax;

void test_error(void);
void test_math_lsqrt(void);
void test_vecl(void);
void test_vecld(void);
void test_vecll(void);
void test_matld(void);
void test_matld_vecld(void);
void test_matll(void);
void test_matld_lowrank(void);
void test_matld_drop(void);
void test_matll_drop(void);
void test_matld_glq(void);
void test_matld_ipstab(void);
void test_matll_ipstab(void);
void test_math_gaussll(void);
void test_math_legendrell(void);
void test_fmmld(void);
void test_fmmld_preproc(fxt_fmmld *fmm);
void test_fmmld_estimate(fxt_fmmld *fmm);
void test_fmmld_regions(fxt_fmmld *fmm);
void test_fmmld_network(fxt_fmmld *fmm);
void test_fmmld_matrix(fxt_fmmld *fmm);
long test_fxtld(long p, long n, long m, double eps);
void test_flptld(char*);
void test_faltld(char*);
void test_faltld_comp(fxt_faltld*);

#endif
