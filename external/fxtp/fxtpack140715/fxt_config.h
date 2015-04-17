#ifndef __FXT_CONFIG_H_INCLUDED__
#define __FXT_CONFIG_H_INCLUDED__

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  parametric choices for FXTPACK
*********************************************************************/

/* The preprocessed data file is also used for FORTRAN routines,
 so the format is that of 'unformatted' file of FORTRAN.  To read
 the file from FORTRAN, the sizes of integers and pointers must
 coincide.  INTEGER_TYPE and DATSIZE_TYPE specifies the type of
 integers and pointers in C.  Values are 1 for short, 2 for int,
 and 3 for long. */

#define INTEGER_TYPE 2
#define DATSIZE_TYPE 3

/* The precision of the resulting transform may not attend the
 required precision.  Major reasons are rounding errors and the
 convergence errors of the power method.  Set the following macro
 zero to prevent aborting preprocessing even when larger error
 than the required precision is found. */

#define DISALLOWHIGHERROR 1

/* Set the following macro non-zero to print progress of the
 preprocessing.  Preprocessing may takes very long time, so it
 may help you to check whether it actually runs or not. */

#define SHOWFXTMAKE 0

/* FXTPACK evaluates errors in 2-norm by using power method. Since
 it is iterative method, one must define convergence criteria.
 Requesting too high precision causes almost useless increase of
 computational time, but too low precision may lead to unstable
 results. */

#define POWERPREC 0.01

/* Next determines the least number of power method iteration */

#define MINPOWERITER 10

/* Power method may fail to converge to the largest eigenvalue,
 so FXTPACK executes power method several times until the maximum
 eigenvalue found twice.  However this scheme fails when the value
 is near the round-off error level.  So FXTPACK stops repeating
 power method after MAXPOWRTRY trials. */

#define MAXPOWERTRY 10

#endif
