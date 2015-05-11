#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_config.h"
#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  file I/O for FXTPACK
*********************************************************************/

#if INTEGER_TYPE == 1
typedef short integer;
#elif INTEGER_TYPE == 2
typedef int integer;
#else
typedef long integer;
#endif

#if DATSIZE_TYPE == 1
typedef short datsize;
#elif DATSIZE_TYPE == 2
typedef int datsize;
#else
typedef long datsize;
#endif

void fxt_file_writestr(FILE *f, const char *s) {
  size_t nw;  long n;  datsize dsz;

  /* check null pointers */
  if (f == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_writestr: null pointer\n");
    return;
  }

  n = strlen(s) + 1;
  dsz = n * sizeof(char);

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writestr: error in writing header\n");
    return;
  }

  nw = fwrite(s, sizeof(char), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writestr: error in writing data\n");
    return;
  }

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writestr: error in writing tailer\n");
    return;
  }
}


void fxt_file_readstr(FILE *f, char *s, long n) {
  size_t nw;  datsize dsz, dsz2;

  /* check null pointers */
  if (f == NULL || s == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readstr: null pointer\n");
    return;
  }

  nw = fread(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readstr: error in reading header\n");
    return;
  }

  if (dsz != n * sizeof(char)) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readstr: data count mismatch\n");
    return;
  }

  nw = fread(s, sizeof(char), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readstr: error in reading data\n");
    return;
  }

  nw = fread(&dsz2, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readstr: error in reading tailer\n");
    return;
  }

  if (dsz != dsz2) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readstr: header/tailer mismatch\n");
    return;
  }
}


void fxt_file_writelongs(FILE *f, const long *x, long n) {
  size_t nw;  datsize dsz;
  integer *xx;  long i;

  /* check null pointers */
  if (f == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_writelongs: null pointer\n");
    return;
  }

  dsz = n * sizeof(integer);

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writelongs: error in writing header\n");
    return;
  }

  /* convert to integer */
  xx = (integer*) malloc(sizeof(integer) * n);
  if (xx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writelongs: allocation failed\n");
    return;
  }

  for (i=0; i< n; i++)
    xx[i] = (integer) x[i];

  nw = fwrite(xx, sizeof(integer), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writelongs: error in writing data\n");
    return;
  }

  free(xx);

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writelongs: error in writing tailer\n");
    return;
  }
}


void fxt_file_readlongs(FILE *f, long *x, long n) {
  size_t nw;  datsize dsz, dsz2;
  integer *xx;  long i;

  /* check null pointers */
  if (f == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readlongs: null pointer\n");
    return;
  }

  nw = fread(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readlongs: error in reading header\n");
    return;
  }

  if (dsz != n * sizeof(integer)) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readlongs: data count mismatch\n");
    return;
  }

  xx = (integer*) malloc(sizeof(integer) * n);
  if (xx == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readlongs: allocation failed\n");
    return;
  }

  nw = fread(xx, sizeof(integer), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readlongs: error in reading data\n");
    return;
  }

  for (i=0; i< n; i++)
    x[i] = (long) xx[i];

  free(xx);

  nw = fread(&dsz2, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readlongs: error in reading tailer\n");
    return;
  }
}


void fxt_file_writedoubles(FILE *f, const double *x, long n) {
  size_t nw;  datsize dsz;

  /* check null pointers */
  if (f == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_writedoubles: null pointer\n");
    return;
  }

  dsz = n * sizeof(double);

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writedoubles: error in writing header\n");
    return;
  }

  nw = fwrite(x, sizeof(double), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writedoubles: error in writing data\n");
    return;
  }

  nw = fwrite(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_writedoubles: error in writing tailer\n");
    return;
  }
}


void fxt_file_readdoubles(FILE *f, double *x, long n) {
  size_t nw;  datsize dsz, dsz2;

  /* check null pointers */
  if (f == NULL || x == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readdoubles: null pointer\n");
    return;
  }

  nw = fread(&dsz, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readdoubles: error in reading header\n");
    return;
  }

  if (dsz != n * sizeof(double)) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_file_readdoubles: data count mismatch\n");
    return;
  }

  nw = fread(x, sizeof(double), n, f);
  if (nw != n) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readdoubles: error in reading data\n");
    return;
  }

  nw = fread(&dsz2, sizeof(datsize), 1, f);
  if (nw != 1) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_file_readdoubles: error in reading tailer\n");
    return;
  }
}
