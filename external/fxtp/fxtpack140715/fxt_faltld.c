#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

#include "fxt_error.h"
#include "fxt_faltld.h"
#include "fxt_faltld_loc.h"

#include "fxt_file.h"

/*********************************************************************
  FXTPACK
   by Reiji Suda, the University of Tokyo
   2014.07.15

  load and delete Fast Associated Legendre Transform
*********************************************************************/

const char* fxt_faltld_header = "FXTPACK FALTLD ver20140715";

/*** load fast associated Legendre transform ***/
fxt_faltld* fxt_faltld_load(char *fname) {
  fxt_faltld *falt;  char s[128];
  long size[3], nm, i;
  FILE *fin;

  /* check pointer */
  if (fname == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_load: null pointer\n");
    return NULL;
  }

  /* open file */
  fin = fopen(fname, "rb");
  if (fin == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_load: cannot open %s\n", fname);
    return NULL;
  }

  /* check header (version) */
  fxt_file_readstr(fin, s, strlen(fxt_faltld_header) + 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (strcmp(s, fxt_faltld_header) != 0) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_load: header (version?) mismatch\n");
    return NULL;
  }

  /* allocate body */
  falt = (fxt_faltld*) malloc(sizeof(fxt_faltld));
  if (falt == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_faltld_load: allocation failed\n");
    return NULL;
  }

  /* read sizes */
  fxt_file_readlongs(fin, size, 3);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;
  
  falt->p = size[0];
  falt->n = size[1];
  nm = size[2];

  /* read precision */
  fxt_file_readdoubles(fin, &falt->prec, 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  /* load Gaussian points */
  falt->x = fxt_vecld_restore(fin);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (falt->x->n != falt->p) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_faltld_load: size mismatch\n");
    return NULL;
  }

  /* load Gaussian weights */
  falt->w = fxt_vecld_restore(fin);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (falt->w->n != falt->p) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_faltld_load: size mismatch\n");
    return NULL;
  }

  /* allocate array of Legendre transforms */
  falt->alt = (faltld_alt**)
    malloc(sizeof(faltld_alt*) * (falt->n + 1));
  if (falt->alt == NULL) {
    fxt_error_set(FXT_ERROR_SYSTEM,
		  "fxt_faltld_load: allocation failed\n");
    return NULL;
  }

  /* allocate order array */
  falt->mv = fxt_vecl_new(nm);

  /* initialize Legendre transform pointers */
  for (i=0; i<= falt->n; i++)
    falt->alt[i] = NULL;

  /* load Legendre transforms */
  for (i=0; i< nm; i++) {
    faltld_alt* alt;

    /* allocate new Legendre transform */
    alt = (faltld_alt*) malloc(sizeof(faltld_alt));
    if (alt == NULL) {
      fxt_error_set(FXT_ERROR_SYSTEM,
		    "fxt_faltld_load: allocation failed\n");
      return NULL;
    }

    /* read order m */
    fxt_file_readlongs(fin, &alt->m, 1);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    falt->mv->v[i] = alt->m;

    /* check m */
    if (alt->m < 0 || falt->n < alt->m) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_faltld_load: irregal order m\n");
      return NULL;
    }

    if (falt->alt[alt->m] != NULL) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_faltld_load: duplicated order m\n");
      return NULL;
    }

    /* load even part */
    alt->ear = fxt_sarld_restore(fin);
    if (fxt_error_raise() > FXT_ERROR_WARN)
      return NULL;

    if (alt->ear->ncol != (falt->n - alt->m + 2) / 2 ||
	alt->ear->nrow != (falt->p + 1) / 2) {
      fxt_error_set(FXT_ERROR_FXTBUG,
		    "fxt_faltld_load: sizemismatch\n");
      return NULL;
    }

    /* load odd part */
    if (falt->n != alt->m) {
      alt->oar = fxt_sarld_restore(fin);
      if (fxt_error_raise() > FXT_ERROR_WARN)
	return NULL;

      if (alt->oar->ncol != (falt->n - alt->m + 1) / 2 ||
	  alt->oar->nrow != falt->p / 2) {
	fxt_error_set(FXT_ERROR_FXTBUG,
		      "fxt_faltld_load: size mismatch\n");
	return NULL;
      }
    } else
      alt->oar = NULL;

    /* connect to the body */
    alt->par = falt;

    falt->alt[alt->m] = alt;
  }

  /* check tailer */
  fxt_file_readstr(fin, s, strlen(fxt_faltld_header) + 1);
  if (fxt_error_raise() > FXT_ERROR_WARN)
    return NULL;

  if (strcmp(s, fxt_faltld_header) != 0) {
    fxt_error_set(FXT_ERROR_FXTBUG,
		  "fxt_faltld_load: tailer mismatch\n");
    return NULL;
  }

  /* close file */
  fclose(fin);

  return falt;
}



void fxt_faltld_del(fxt_faltld *falt) {
  long i;

  /* check pointers */
  if (falt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_del: null pointer\n");
    return;
  }

  for (i= falt->mv->n - 1; i >= 0; i --) {

    /* the data to be deallocated */
    faltld_alt *alt = falt->alt[falt->mv->v[i]];

    /* deallocate odd sarld */
    if (alt->oar != NULL) {
      fxt_sarld_del(alt->oar);
      alt->oar = NULL;
    }

    /* deallocate even sarld */
    fxt_sarld_del(alt->ear);
    alt->ear = NULL;

    /* nullify pointers */
    alt->par = NULL;

    /* deallocate it */
    free(alt);
  }

  /* deallocate array of order m */
  fxt_vecl_del(falt->mv);
  falt->mv = NULL;

  /* deallocate array of Legendre transforms */
  free(falt->alt);
  falt->alt = NULL;

  /* deallocate Gaussian weights */
  fxt_vecld_del(falt->w);
  falt->w = NULL;

  /* deallocate Gaussian nodes */
  fxt_vecld_del(falt->x);
  falt->x = NULL;

  /* deallocate the body */
  free(falt);
}


/*** required working array size ***/
long fxt_faltld_wsize(fxt_faltld *falt, long m) {
  long wsize;
  fxt_sarld *ar;

  /* check pointers */
  if (falt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_wsize: null pointer\n");
    return 0;
  }

  /* check order m */
  if (m < 0 || falt->n < m || falt->alt[m] == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_wsize: irregal order m\n");
    return 0;
  }

  /* even part */
  ar = falt->alt[m]->ear;

  wsize = ar->ncol + ar->nrow + ar->ntmp;

  /* odd part */
  if (m == falt->n)
    return wsize;

  ar = falt->alt[m]->oar;

  if (wsize < ar->ncol + ar->nrow + ar->ntmp)
    wsize = ar->ncol + ar->nrow + ar->ntmp;

  return wsize;
}


/*** maximum workin array size ***/
long fxt_faltld_wsizemax(fxt_faltld *falt) {
  long i, wsize;

  /* check pointer */
  if (falt == NULL) {
    fxt_error_set(FXT_ERROR_USAGE,
		  "fxt_faltld_wsizemax: null pointer\n");
    return 0;
  }

  wsize = 0;
  for (i=0; i< falt->mv->n; i++) {
    long ws = fxt_faltld_wsize(falt, falt->mv->v[i]);

    if (fxt_error_raise() > FXT_ERROR_WARN)
      return 0;

    if (ws > wsize)
      wsize = ws;
  }

  return wsize;
}
