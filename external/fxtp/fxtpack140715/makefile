#  FXTPACK
#   by Reiji Suda, the University of Tokyo
#   2014.07.15

default : fxtpp_faltld sample_faltld fxtpp_flptld sample_flptld\
 test_fxtpack

# --- GNU C Compiler
CC = gcc -Wall -O4 -pedantic -ansi -Wundef

LOBJ = fxt_error.o fxt_file.o fxt_vecl.o fxt_vecld.o\
 fxt_sarld.o fxt_sarld_comp.o\
 fxt_faltld.o fxt_faltld_comp.o fxt_flptld.o fxt_flptld_comp.o

OBJ = fxt_error.o fxt_file.o\
 fxt_math_lsqrt.o fxt_math_gaussll.o fxt_math_legendrell.o\
 fxt_vecl.o fxt_vecld.o fxt_vecll.o fxt_matld.o fxt_matll.o\
 fxt_matld_vecld.o fxt_matld_lowrank.o fxt_matld_drop.o\
 fxt_matld_glq.o fxt_matld_ipstab.o fxt_matld_lagld.o\
 fxt_matll_drop.o fxt_matll_ipstab.o\
 fxt_fmmld.o fxt_fmmld_preproc.o\
 fxt_fmmld_regions.o fxt_fmmld_network.o fxt_fmmld_matrix.o\
 fxt_fmmld_estimate.o fxt_fmmld_evaluate.o fxt_fmmld_lagld.o\
 fxt_fmmld_mulexpmat.o fxt_fmmld_locexpmat.o fxt_fmmld_scalevec.o\
 fxt_fmmld_evl.o fxt_fmmld_exp.o fxt_fmmld_scalefmm.o\
 fxt_fxtld.o fxt_fxtld_make.o fxt_fxtld_fmm.o fxt_fxtld_div.o\
 fxt_fxtld_amat.o fxt_fxtld_info.o fxt_fxtld_comp.o\
 fxt_fxtld_scale.o fxt_fxtld_lagld.o\
 fxt_lagld.o fxt_lagld_oper.o fxt_lagld_comp.o\
 fxt_lagld_drop.o fxt_lagld_normal.o\
 fxt_sarld.o fxt_sarld_lagld.o\
 fxt_sarld_comp.o fxt_sarld_error.o fxt_sarld_scale.o\
 fxt_flptld.o fxt_flptld_preproc.o fxt_flptld_comp.o\
 fxt_faltld.o fxt_faltld_preproc.o fxt_faltld_comp.o

fxtpp_faltld : fxtpp_faltld.o $(OBJ)
	$(CC) -o fxtpp_faltld fxtpp_faltld.o $(OBJ) -lm

sample_faltld : sample_faltld.o $(LOBJ)
	$(CC) -o sample_faltld sample_faltld.o $(LOBJ) -lm

fxtpp_flptld : fxtpp_flptld.o $(OBJ)
	$(CC) -o fxtpp_flptld fxtpp_flptld.o $(OBJ) -lm

sample_flptld : sample_flptld.o $(LOBJ)
	$(CC) -o sample_flptld sample_flptld.o $(LOBJ) -lm
fxt_error.o: fxt_error.c fxt_error.h
fxt_faltld.o: fxt_faltld.c fxt_error.h fxt_faltld.h \
  fxt_vecl.h fxt_vecld.h fxt_faltld_loc.h fxt_sarld.h fxt_lagld.h \
  fxt_lagld_loc.h fxt_matld.h fxt_file.h
fxt_faltld_comp.o: fxt_faltld_comp.c fxt_error.h \
  fxt_faltld.h fxt_vecl.h fxt_vecld.h fxt_faltld_loc.h fxt_sarld.h \
  fxt_lagld.h fxt_lagld_loc.h fxt_matld.h
fxt_faltld_preproc.o: fxt_faltld_preproc.c fxt_error.h \
  fxt_config.h fxt_faltld.h fxt_vecl.h fxt_vecld.h fxt_faltld_loc.h \
  fxt_sarld.h fxt_lagld.h fxt_lagld_loc.h fxt_matld.h \
  fxt_math.h fxt_vecll.h fxt_matll.h fxt_fxtld.h fxt_file.h
fxt_file.o: fxt_file.c fxt_error.h fxt_config.h fxt_file.h
fxt_flptld.o: fxt_flptld.c fxt_error.h fxt_flptld.h \
  fxt_vecld.h fxt_vecl.h fxt_sarld.h fxt_lagld.h fxt_lagld_loc.h \
  fxt_matld.h fxt_flptld_loc.h fxt_file.h
fxt_flptld_comp.o: fxt_flptld_comp.c fxt_error.h \
  fxt_flptld.h fxt_vecld.h fxt_vecl.h fxt_sarld.h fxt_lagld.h \
  fxt_lagld_loc.h fxt_matld.h fxt_flptld_loc.h
fxt_flptld_preproc.o: fxt_flptld_preproc.c fxt_error.h \
  fxt_config.h fxt_flptld.h fxt_vecld.h fxt_vecl.h fxt_sarld.h \
  fxt_lagld.h fxt_lagld_loc.h fxt_matld.h fxt_flptld_loc.h fxt_math.h \
  fxt_vecll.h fxt_matll.h fxt_fxtld.h fxt_file.h
fxt_fmmld.o: fxt_fmmld.c fxt_error.h fxt_config.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_estimate.o: fxt_fmmld_estimate.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_evaluate.o: fxt_fmmld_evaluate.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_evl.o: fxt_fmmld_evl.c fxt_error.h fxt_fmmld.h \
  fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_fmmld_loc.h \
  fxt_lagld_loc.h
fxt_fmmld_exp.o: fxt_fmmld_exp.c fxt_error.h fxt_fmmld.h \
  fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_fmmld_loc.h \
  fxt_lagld_loc.h
fxt_fmmld_lagld.o: fxt_fmmld_lagld.c fxt_error.h \
  fxt_fmmld_loc.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_lagld_loc.h
fxt_fmmld_locexpmat.o: fxt_fmmld_locexpmat.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_matrix.o: fxt_fmmld_matrix.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_mulexpmat.o: fxt_fmmld_mulexpmat.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_network.o: fxt_fmmld_network.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_preproc.o: fxt_fmmld_preproc.c fxt_error.h \
  fxt_config.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_regions.o: fxt_fmmld_regions.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_scalefmm.o: fxt_fmmld_scalefmm.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fmmld_scalevec.o: fxt_fmmld_scalevec.c fxt_error.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h \
  fxt_fmmld_loc.h fxt_lagld_loc.h
fxt_fxtld.o: fxt_fxtld.c fxt_error.h fxt_config.h \
  fxt_fxtld.h fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h \
  fxt_lagld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_amat.o: fxt_fxtld_amat.c fxt_error.h fxt_config.h \
  fxt_fxtld.h fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h \
  fxt_lagld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_comp.o: fxt_fxtld_comp.c fxt_error.h fxt_fxtld.h \
  fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h \
  fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_div.o: fxt_fxtld_div.c fxt_error.h fxt_fxtld.h \
  fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h \
  fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_fmm.o: fxt_fxtld_fmm.c fxt_error.h fxt_config.h \
  fxt_fxtld.h fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h \
  fxt_lagld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_info.o: fxt_fxtld_info.c fxt_error.h fxt_config.h \
  fxt_fxtld.h fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h \
  fxt_lagld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_lagld.o: fxt_fxtld_lagld.c fxt_error.h fxt_lagld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld_loc.h fxt_fxtld.h fxt_matll.h \
  fxt_vecll.h fxt_matld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_make.o: fxt_fxtld_make.c fxt_error.h fxt_config.h \
  fxt_fxtld.h fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h \
  fxt_lagld.h fxt_fxtld_loc.h fxt_fmmld.h
fxt_fxtld_scale.o: fxt_fxtld_scale.c fxt_error.h fxt_fxtld.h \
  fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h \
  fxt_fxtld_loc.h fxt_fmmld.h
fxt_lagld.o: fxt_lagld.c fxt_error.h fxt_lagld.h fxt_vecld.h \
  fxt_vecl.h fxt_lagld_loc.h
fxt_lagld_comp.o: fxt_lagld_comp.c fxt_error.h fxt_lagld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld_loc.h
fxt_lagld_drop.o: fxt_lagld_drop.c fxt_error.h fxt_lagld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld_loc.h
fxt_lagld_normal.o: fxt_lagld_normal.c fxt_error.h \
  fxt_lagld.h fxt_vecld.h fxt_vecl.h fxt_lagld_loc.h
fxt_lagld_oper.o: fxt_lagld_oper.c fxt_error.h fxt_lagld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld_loc.h
fxt_math_gaussll.o: fxt_math_gaussll.c fxt_error.h \
  fxt_math.h fxt_vecll.h fxt_vecl.h fxt_vecld.h fxt_matll.h fxt_matld.h \
  fxt_lagld.h
fxt_math_legendrell.o: fxt_math_legendrell.c fxt_error.h \
  fxt_math.h fxt_vecll.h fxt_vecl.h fxt_vecld.h fxt_matll.h fxt_matld.h \
  fxt_lagld.h
fxt_math_lsqrt.o: fxt_math_lsqrt.c fxt_error.h fxt_math.h \
  fxt_vecll.h fxt_vecl.h fxt_vecld.h fxt_matll.h fxt_matld.h fxt_lagld.h
fxt_matld.o: fxt_matld.c fxt_error.h fxt_config.h \
  fxt_matld.h fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matld_drop.o: fxt_matld_drop.c fxt_error.h fxt_matld.h \
  fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matld_glq.o: fxt_matld_glq.c fxt_error.h fxt_matld.h \
  fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matld_ipstab.o: fxt_matld_ipstab.c fxt_error.h \
  fxt_matld.h fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matld_lagld.o: fxt_matld_lagld.c fxt_error.h \
  fxt_config.h fxt_matld.h fxt_vecl.h fxt_vecld.h fxt_lagld.h \
  fxt_lagld_loc.h
fxt_matld_lowrank.o: fxt_matld_lowrank.c fxt_error.h \
  fxt_matld.h fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matld_vecld.o: fxt_matld_vecld.c fxt_error.h fxt_matld.h \
  fxt_vecl.h fxt_vecld.h fxt_lagld.h
fxt_matll.o: fxt_matll.c fxt_error.h fxt_matll.h fxt_vecl.h \
  fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h
fxt_matll_drop.o: fxt_matll_drop.c fxt_error.h fxt_matll.h \
  fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h
fxt_matll_ipstab.o: fxt_matll_ipstab.c fxt_error.h \
  fxt_matll.h fxt_vecl.h fxt_vecll.h fxt_vecld.h fxt_matld.h fxt_lagld.h
fxt_sarld.o: fxt_sarld.c fxt_error.h fxt_sarld.h fxt_vecld.h \
  fxt_vecl.h fxt_lagld.h fxt_lagld_loc.h fxt_matld.h fxt_file.h
fxt_sarld_comp.o: fxt_sarld_comp.c fxt_error.h fxt_sarld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld.h fxt_lagld_loc.h fxt_matld.h
fxt_sarld_error.o: fxt_sarld_error.c fxt_error.h \
  fxt_config.h fxt_sarld.h fxt_vecld.h fxt_vecl.h fxt_lagld.h \
  fxt_lagld_loc.h fxt_matld.h
fxt_sarld_lagld.o: fxt_sarld_lagld.c fxt_error.h fxt_sarld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld.h fxt_lagld_loc.h fxt_matld.h
fxt_sarld_scale.o: fxt_sarld_scale.c fxt_error.h fxt_sarld.h \
  fxt_vecld.h fxt_vecl.h fxt_lagld.h fxt_lagld_loc.h fxt_matld.h
fxt_vecl.o: fxt_vecl.c fxt_error.h fxt_vecl.h
fxt_vecld.o: fxt_vecld.c fxt_error.h fxt_vecld.h fxt_vecl.h \
  fxt_file.h
fxt_vecll.o: fxt_vecll.c fxt_error.h fxt_vecll.h fxt_vecl.h \
  fxt_vecld.h fxt_math.h fxt_matll.h fxt_matld.h fxt_lagld.h
fxtpp_faltld.o: fxtpp_faltld.c fxt_error.h fxt_faltld.h \
  fxt_vecl.h fxt_vecld.h
fxtpp_flptld.o: fxtpp_flptld.c fxt_error.h fxt_flptld.h \
  fxt_vecld.h fxt_vecl.h fxt_sarld.h fxt_lagld.h fxt_lagld_loc.h \
  fxt_matld.h
sample_faltld.o: sample_faltld.c fxt_error.h fxt_faltld.h \
  fxt_vecl.h fxt_vecld.h
sample_flptld.o: sample_flptld.c fxt_error.h fxt_flptld.h \
  fxt_vecld.h fxt_vecl.h fxt_sarld.h fxt_lagld.h fxt_lagld_loc.h \
  fxt_matld.h

TOBJ = test_error.o\
 test_math_lsqrt.o test_math_gaussll.o test_math_legendrell.o\
 test_vecl.o test_vecld.o test_vecll.o test_matld.o test_matll.o\
 test_matld_vecld.o test_matld_drop.o test_matld_lowrank.o\
 test_matld_glq.o test_matld_ipstab.o\
 test_matll_drop.o test_matll_ipstab.o\
 test_fmmld.o test_fmmld_preproc.o test_fmmld_regions.o\
 test_fmmld_network.o test_fmmld_matrix.o test_fmmld_estimate.o\
 test_fxtld.o test_flptld.o test_faltld.o test_faltld_comp.o

test_fxtpack : test_fxtpack.o $(OBJ) $(TOBJ)
	$(CC) -o test_fxtpack test_fxtpack.o $(OBJ) $(TOBJ) -lm
test_error.o: test_error.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h
test_faltld.o: test_faltld.c fxt_error.h fxt_config.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_faltld_loc.h fxt_math.h fxt_vecll.h fxt_matll.h
test_faltld_comp.o: test_faltld_comp.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_flptld.o: test_flptld.c fxt_error.h fxt_config.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_flptld.h fxt_flptld_loc.h fxt_math.h fxt_vecll.h fxt_matll.h
test_fmmld.o: test_fmmld.c fxt_error.h fxt_config.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fmmld_estimate.o: test_fmmld_estimate.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fmmld_matrix.o: test_fmmld_matrix.c fxt_error.h \
  fxt_config.h test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h \
  fxt_matld.h fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fmmld_network.o: test_fmmld_network.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fmmld_preproc.o: test_fmmld_preproc.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fmmld_regions.o: test_fmmld_regions.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fmmld_loc.h
test_fxtld.o: test_fxtld.c fxt_error.h fxt_config.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_fxtld.h fxt_matll.h fxt_vecll.h fxt_fxtld_loc.h fxt_math.h
test_fxtpack.o: test_fxtpack.c test_fxtpack.h fxt_fmmld.h fxt_vecld.h \
  fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h fxt_sarld.h \
  fxt_lagld_loc.h fxt_error.h fxt_flptld.h
test_math_gaussll.o: test_math_gaussll.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_math.h fxt_vecll.h fxt_matll.h
test_math_legendrell.o: test_math_legendrell.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_math.h fxt_vecll.h fxt_matll.h
test_math_lsqrt.o: test_math_lsqrt.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_math.h fxt_vecll.h fxt_matll.h
test_matld.o: test_matld.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h
test_matld_drop.o: test_matld_drop.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_matld_glq.o: test_matld_glq.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_matld_ipstab.o: test_matld_ipstab.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_matld_lowrank.o: test_matld_lowrank.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_matld_vecld.o: test_matld_vecld.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h
test_matll.o: test_matll.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h fxt_matll.h fxt_vecll.h
test_matll_drop.o: test_matll_drop.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_matll.h fxt_vecll.h
test_matll_ipstab.o: test_matll_ipstab.c fxt_error.h \
  test_fxtpack.h fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h \
  fxt_lagld.h fxt_faltld.h fxt_sarld.h fxt_lagld_loc.h \
  fxt_matll.h fxt_vecll.h
test_vecl.o: test_vecl.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h
test_vecld.o: test_vecld.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h
test_vecll.o: test_vecll.c fxt_error.h test_fxtpack.h \
  fxt_fmmld.h fxt_vecld.h fxt_vecl.h fxt_matld.h fxt_lagld.h fxt_faltld.h \
  fxt_sarld.h fxt_lagld_loc.h fxt_vecll.h
