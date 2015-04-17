#ifndef __FXT_FILE_H_INCLUDED__
#define __FXT_FILE_H_INCLUDED__

#include <stdio.h>

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  file I/O routines for FXTPACK
*********************************************************************/

/*** write a string ***/
void fxt_file_writestr(FILE *f, const char *s);

/*** read a string ***/
void fxt_file_readstr(FILE *f, char *s, long size);

/*** write an array of long ***/
void fxt_file_writelongs(FILE *f, const long *v, long n);

/*** read an array of long ***/
void fxt_file_readlongs(FILE *f, long *v, long n);

/*** write an array of double ***/
void fxt_file_writedoubles(FILE *f, const double *v, long n);

/*** read an array of double ***/
void fxt_file_readdoubles(FILE *f, double *v, long n);

#endif
