#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_flptld.h"
#include "fxt_flptld_loc.h"

#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  load/deallocate Fast Legendre Polynomial Transform data stucture
*********************************************************************/

const char* fxt_flptld_header = "FXTPACK FLPTLD ver20140715";

/*** load fast Legendre polynomial transform ***/
fxt_flptld* fxt_flptld_load(char *fname) {
  fxt_flptld *flpt;  char s[128];
  long size[2];  FILE *fin;

  /* check null pointers */
  if (fname == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_load: null pointer\n");
    return NULL;
  }

  /* open file */
  fin = fopen(fname, "rb");
  if (fin == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_load: cannot open %s\n", fname);
    return NULL;
  }

  /* check header (version) */
  fxt_file_readstr(fin, s, strlen(fxt_flptld_header) + 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (strcmp(s, fxt_flptld_header) != 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_load: header (version?) mismatch\n");
    return NULL;
  }

  /* allocate body */
  flpt = (fxt_flptld*) malloc(sizeof(fxt_flptld));
  if (flpt == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_flptld_load: allocation failed\n");
    return NULL;
  }

  /* read sizes */
  fxt_file_readlongs(fin, size, 2);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;
  
  flpt->p = size[0];
  flpt->n = size[1];

  /* read precision */
  fxt_file_readdoubles(fin, &flpt->prec, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* load Gaussian points */
  flpt->x = fxt_vecld_restore(fin);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (flpt->x->n != flpt->p) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_flptld_load: size mismatch\n");
    return NULL;
  }

  /* load Gaussian weights */
  flpt->w = fxt_vecld_restore(fin);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (flpt->w->n != flpt->p) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_flptld_load: size mismatch\n");
    return NULL;
  }

  /* load even part */
  flpt->ear = fxt_sarld_restore(fin);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (flpt->ear->ncol != (flpt->n + 2) / 2 ||
      flpt->ear->nrow != (flpt->p + 1) / 2) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_flptld_load: sizemismatch\n");
    return NULL;
  }

  /* load odd part */
  if (flpt->n != 0) {
    flpt->oar = fxt_sarld_restore(fin);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    if (flpt->oar->ncol != (flpt->n + 1) / 2 ||
	flpt->oar->nrow != flpt->p / 2) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_flptld_load: size mismatch\n");
      return NULL;
    }
  } else
    flpt->oar = NULL;

  /* check tailer */
  fxt_file_readstr(fin, s, strlen(fxt_flptld_header) + 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (strcmp(s, fxt_flptld_header) != 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_flptld_load: tailer mismatch\n");
    return NULL;
  }

  /* close file */
  fclose(fin);

  return flpt;
}


void fxt_flptld_del(fxt_flptld *flpt) {

  /* check null pointers */
  if (flpt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_del: null pointer\n");
    return;
  }

  if (flpt->oar != NULL) {
    fxt_sarld_del(flpt->oar);
    flpt->oar = NULL;
  }

  fxt_sarld_del(flpt->ear);
  flpt->ear = NULL;

  fxt_vecld_del(flpt->w);
  flpt->w = NULL;

  fxt_vecld_del(flpt->x);
  flpt->x = NULL;

  free(flpt);
}


/*** compute the working array size ***/
long fxt_flptld_wsize(fxt_flptld *flpt) {
  long wsize;
  fxt_sarld *ar;

  /* check pointers */
  if (flpt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_flptld_wsize: null pointer\n");
    return 0;
  }

  /* even part */
  ar = flpt->ear;

  wsize = ar->ncol + ar->nrow + ar->ntmp;

  /* odd part */
  if (flpt->n == 0)
    return wsize;

  ar = flpt->oar;

  if (wsize < ar->ncol + ar->nrow + ar->ntmp)
    wsize = ar->ncol + ar->nrow + ar->ntmp;

  return wsize;
}
